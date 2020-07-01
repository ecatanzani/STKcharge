import sys
from datetime import datetime
from argparse import ArgumentParser


start = datetime.now()


def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Acceptance Adder")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input condor jobs WD')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output ROOT TFile')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')
    parser.add_argument("-c", "--check", dest='check', default=False,
                        action='store_true', help='check tmp acceptance ROOT files')
    
    opts = parser.parse_args(args)

    # Load analysis functions
    sys.path.append("moduls")
    from adder import compute_final_histos
    from scanDirs import getListOfFiles
    from submitJobs import resubmit_condor_jobs

    good_dirs, skipped_dirs, skipped_file_notFinalDir, skipped_file_notROOTfile, skipped_file_notReadable, skipped_file_noKeys= getListOfFiles(opts.input)

    if opts.verbose:
        print('Found {} GOOD condor directories'.format(len(good_dirs)))

    if skipped_dirs:
        print('Found {} BAD condor directories...\n'.format(len(skipped_dirs)))
        print('Found {} directories with no output folder'.format(skipped_file_notFinalDir))
        print('Found {} directories with no output ROOT file'.format(skipped_file_notROOTfile))
        print('Found {} directories with corrupted output ROOT file'.format(skipped_file_notReadable))
        print('Found {} directories where output ROOT file has no keys\n'.format(skipped_file_noKeys))
        print('Here the folders list...\n')

        for idx, elm in enumerate(skipped_dirs):
            print('Skipped {} directory: {}'.format(idx, elm))

        print('\nResubmitting HTCondor jobs for {} directories\n'.format(len(skipped_dirs)))
        resubmit_condor_jobs(skipped_dirs, opts)

    else:
        if opts.verbose:
            print('Found {} BAD condor directories...'.format(len(skipped_dirs)))

    if not opts.check:
        if opts.output:
            compute_final_histos(good_dirs, opts)
        else:
            print("!!! Missing output ROOT file... please specify a path")

if __name__ == "__main__":
    main()
    print(datetime.now()-start)
