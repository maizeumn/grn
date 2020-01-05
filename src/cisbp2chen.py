#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import os.path as op

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'convert cis_bp PWM file to chen format'
    )
    parser.add_argument('fi', help = 'input cis_bp PWM file')
    parser.add_argument('--id', default = '', help = 'motif id [file basename]')
    args = parser.parse_args()

    mid = args.id
    if mid == '':
        mid = op.basename(op.splitext(args.fi)[0])

    fhi = open(args.fi, 'r')
    for line in fhi:
        ps = line.strip().split("\t")
        if len(ps) != 5:
            print("not 5 fields: %s" % line)
            sys.exit(1)
        if ps[0].lower() == 'pos':
            print(">%s" % mid)
            continue
        print("\t".join(ps[1:5]))
    fhi.close()
