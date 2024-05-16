#!/usr/bin/env python

import argparse

import parse_protocol
import pyfastx
import utils

logger = utils.get_logger(__name__)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--whitelist")
    parser.add_argument("--pattern")
    args = parser.parse_args()

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    # protocol
    if args.protocol == "auto":
        runner = parse_protocol.Auto(fq1_list, args.sample, args.assets_dir)
        protocol, protocol_meta = runner.run()
    else:
        protocol = args.protocol
        if protocol == "new":
            protocol_meta = {}
            protocol_meta["bc"] = args.whitelist.split(",")
            protocol_meta["pattern_dict"] = parse_protocol.parse_pattern(args.pattern)
        else:
            protocol_meta = parse_protocol.get_protocol_dict(args.assets_dir)[protocol]

    pattern_dict = protocol_meta["pattern_dict"]
    raw_list, mismatch_list = parse_protocol.get_raw_mismatch(protocol_meta["bc"], 1)

    # out_fq
    out_fq_fn = {x: f"{args.sample}_R{x}.fq.gz" for x in [1, 2]}
    outdict = {k: utils.openfile(v, "wt") for k, v in out_fq_fn.items()}

    n = 0
    for fq1, fq2 in zip(fq1_list, fq2_list):
        logger.info(f"running {fq1} and {fq2}")
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)

        for (name1, seq1, qual1), (name2, seq2, qual2) in zip(fq1, fq2):
            n += 1
            bc_list = [seq2[x] for x in pattern_dict["C"]]
            valid, corrected, corrected_seq = parse_protocol.check_seq_mismatch(bc_list, raw_list, mismatch_list)
            if valid:
                read_name = f"{corrected_seq}:{n}"
                outdict[1].write(utils.fastq_str(read_name, seq1, qual1))
                outdict[2].write(utils.fastq_str(read_name, seq2, qual2))

    # output protocol json
    outf = f"{args.sample}.protocol.json"
    data = {"protocol": protocol}
    utils.write_json(data, outf)
