workflow run_wf {

  take: in_

  main:

    out_ = in_
      | map{ id, f -> [ id.replaceAll("__R[12]", ""), [ input: f.input ] , id ] }
      | groupTuple(sort: "hash")
      | map{ new_id, inputs, ids ->
        r1 = inputs.collect{it.input}.findAll{it =~ "_R1_"}[0]
        r2 = inputs.collect{it.input}.findAll{it =~ "_R2_"}[0]
        [
          new_id,
          [
            r1: r1,
            r2: r2,
            _meta: [ join_id: ids[0] ]
          ]
        ]
      }

  emit: out_

}

