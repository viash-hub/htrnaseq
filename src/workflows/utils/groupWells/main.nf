workflow run_wf {

  take: in_

  main:

    out_ = in_
      | map{ id, state -> [state.pool, state, id]}
      | groupTuple(sort: "hash")
      | map{ new_id, inputs, original_ids ->
        [
          new_id,
          [
            output_r1: inputs.collect{it.input_r1},
            output_r2: inputs.collect{it.input_r2},
            wells: inputs.collect{it.well},
            _meta: [ join_id: original_ids[0] ]
          ]
        ]
      }

  emit: out_

}

