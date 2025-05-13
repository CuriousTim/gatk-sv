version 1.0

import "Structs.wdl"

workflow MakeCohortBatchMap {
  input {
    String sample_set_set_id
    String workspace_namespace
    String workspace_name
    String jq_gcloud_linux_docker
    RuntimeAttr? runtime_override
  }

  call GatherCohortSamples {
    input:
      workspace_namespace = workspace_namespace,
      workspace_name = workspace_name,
      sample_set_set_id = sample_set_set_id,
      jq_gcloud_linux_docker = jq_gcloud_linux_docker,
      runtime_override = runtime_override
  }

  output {
    File batch_map = GatherCohortSamples.batch_map
  }
}

task GatherCohortSamples {
  input {
    String workspace_namespace
    String workspace_name
    String sample_set_set_id
    String jq_gcloud_linux_docker
    RuntimeAttr? runtime_override
  }
  
  RuntimeAttr runtime_default = object {
    cpu_cores: 1,
    mem_gb: 1,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: 16
  }
  RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

  runtime {
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
      cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
      disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
      docker: jq_gcloud_linux_docker
      maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
      memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
      preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
  }

  String sample_batch_map = sample_set_set_id + "-batch_map.tsv"
  command <<<
    set -euo pipefail

    curl --oauth2-bearer "$(gcloud auth application-default print-access-token)" \
      --request GET --header 'accept: */*' \
      'https://api.firecloud.org/api/workspaces/~{workspace_namespace}/~{workspace_name}/entities/sample_set_set/~{sample_set_set_id}' \
      | jq --raw-output '.attributes.sample_sets.items[] | .entityName' > sample_sets.list
    if [[ ! -s sample_sets.list ]]; then
      printf 'no sample sets found\n' >&2
      exit 1
    fi

    read -r id < sample_sets.list
    if [[ "${id}" = 'null' ]]; then
      printf 'no sample sets found\n' >&2
      exit 1
    fi

    declare -i i=1
    : > samples.tsv
    while true; do
      curl --oauth2-bearer "$(gcloud auth application-default print-access-token)" \
        --request GET --header 'accept: application/json' \
        --url-query "page=${i}" \
        --url-query "pageSize=10" \
        --url-query "fields=samples" \
        'https://api.firecloud.org/api/workspaces/~{workspace_namespace}/~{workspace_name}/entityQuery/sample_set' > response.json
      declare -i max_pages
      max_pages=$(jq .resultMetadata.filteredPageCount response.json)
      jq --raw-output \
        '.results[] | .name as $batch | .attributes.samples.items[] | [$batch, .entityName] | @tsv' response.json >> samples.tsv
      if (( i >= max_pages )); then
        break
      fi
      i=$(( i + 1 ))
    done

    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' sample_sets.list samples.tsv > '~{sample_batch_map}'
  >>>

  output {
    File batch_map = sample_batch_map
  }
} 
