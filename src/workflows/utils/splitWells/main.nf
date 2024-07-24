workflow run_wf {

  take: in_

  main:

    out_ = in_

      | flatMap{ id, state ->
        state.input.collect{ p ->
          barcode = (p =~ /.*\\/([ACTG]*|unknown)_R?.*/)[0][1]
          pair_end = (p =~ /.*_(R[12])_.*/)[0][1]
          lane = (id =~ /.*_(L\d+).*/) ? (id =~ /.*_(L\d+).*/)[0][1] : "no_lanes"
          [
            id + "__" + barcode + "__" + pair_end,
            [
              pool: id,
              barcode: barcode,
              barcode_path: p,
              lane: lane,
              pair_end: pair_end,
              _meta: [ join_id: id ]
            ]
          ]
        }
      }

  emit: out_

}

