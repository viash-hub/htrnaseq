workflow run_wf {

  take: in_

  main:

    out_ = in_
      | map{ id, f -> [ id.replaceAll("_L00\\d", "").replaceAll("_lane\\d", ""), [ input_r1: f.input_r1, input_r2: f.input_r2 ] , id ] }
      | groupTuple(sort: "hash")
      | map{ new_id, inputs, ids ->
        [
          new_id,
          [
            output_r1: inputs.collect{it.input_r1},
            output_r2: inputs.collect{it.input_r2},
            _meta: [ join_id: ids[0] ]
          ]
        ]
      }

  emit: out_

}

