#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

import os
import argparse


def main():
    args = parse_arguments()
    convert_links(args.testdir, args.slakodir, args.oldset, args.newset,
                  not args.force)
    if not args.force:
        print('DRY RUN ONLY, CONVERSION WAS NOT DONE')


def parse_arguments():
    msg = 'Converts links slako links from old set to a new one.'
    parser = argparse.ArgumentParser(description=msg)

    msg = 'really do the conversion (by default conversions are only '\
          'shown but not done'
    parser.add_argument('-f', '--force', action='store_true', default=False,
                        help=msg)
    msg = "old set name (e.g. 'mio-0-1')"
    parser.add_argument('oldset', help=msg)
    msg = "new set name (e.g. 'mio-1-1')"
    parser.add_argument('newset', help=msg)
    msg = "directory in which the slako sets are stored (e.g. "\
          "'./autotest/slako')"
    parser.add_argument('slakodir', help=msg)
    msg = "parent directory of the tests (e.g. './autotest')"
    parser.add_argument('testdir', help=msg)
    args = parser.parse_args()
    return args


def convert_links(testdir, slakodir, oldset, newset, dryrun=True):

    slakodir = os.path.realpath(slakodir)
    newslakodir = os.path.join(slakodir, newset)

    for root, dirs, files in os.walk(testdir):
        for ff in files:
            fullpath = os.path.join(root, ff)
            if os.path.islink(fullpath):
                target = os.path.realpath(fullpath)
                setname = os.path.basename(os.path.dirname(target))
                if setname != oldset:
                    continue
                reldir = os.path.relpath(newslakodir, os.path.realpath(root))
                linktarget = os.path.join(reldir, ff)
                print('Relinking:', fullpath, '=>', linktarget)
                if not dryrun:
                    os.remove(fullpath)
                    os.symlink(linktarget, fullpath)



if __name__ == '__main__':
    main()
