import numpy as np
import h5py
import argparse
import os
import time

import surfinBH


def eval_fits(fit, num_tests, function):

    # Gather random params
    q_list = []
    chiA_list = []
    chiB_list = []
    for i in range(num_tests):
        # Generate params randomly within allowed values
        q, chiA, chiB = fit._generate_random_params_for_tests() 
        q_list.append(q)
        chiA_list.append(chiA)
        chiB_list.append(chiB)

    # Test timing
    start_time = time.time()
    for i in range(num_tests):
        # evaluate fit
        y = function(q_list[i], chiA_list[i], chiB_list[i])
    end_time = time.time()

    return (end_time - start_time)/num_tests


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tests evaluation times for'
            ' fits.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", "-n", type=str, required=True, \
            help="Fit name without the surfinBH prefix. Eg. 7dq2.")
    args = parser.parse_args()

    # Load fit
    fit_name = 'surfinBH%s'%args.name
    fit = surfinBH.LoadFits(fit_name)

    # number of test evaluations per fit
    num_tests = 1000

    # remnant mass
    eval_time = eval_fits(fit, num_tests, fit.mf)
    print('Time per evaluation of mf: %.3e s'%eval_time)

    # remnant spin
    eval_time = eval_fits(fit, num_tests, fit.chif)
    print('Time per evaluation of chif: %.3e s'%eval_time)

    # remnant kick
    eval_time = eval_fits(fit, num_tests, fit.vf)
    print('Time per evaluation of vf: %.3e s'%eval_time)
