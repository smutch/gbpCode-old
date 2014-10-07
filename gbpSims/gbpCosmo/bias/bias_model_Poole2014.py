#!/usr/bin/env python

"""Doc-string here..."""

import numpy as np
import argparse

__author__ = "Gregory B. Poole"
__date__ = "2014-05-30"


def set_b_s_Poole_params(V_max, z, halo_type, bias_type):

    """ This function expresses Equations (3), and (4)
        and the values of Table 2 in Poole et al (2014a)

    Args:
        V_max (float) :  halo maximum circular velocity

        z (float) :  redshift of interest

        halo_type (str) : halo type

        bias_type (str) : bias type

    Returns:
        A list containing parameters `b_x`, `b_s` and `V_SF`
    """

    # Set-up a set of dictionaries holding the coefficients of the model
    case = bias_type+'_'+halo_type
    cols = ['real_FoF', 'real_sub', 'zboost_FoF', 'zboost_sub',
            'zspace_FoF', 'zspace_sub']
    V_SF_0_c = [0.02819, 0.05326, 0.31002, 0.31731, 0.20417, 0.21287]
    V_SF_z_c = [-0.13820, -0.16739, -0.20264, -0.15991, -0.29667, -0.22806]
    s_V_0_c = [0.36860, 0.40269, 0.33423, 0.53444, 0.94082, 1.11879]
    s_V_z_c = [0.61547, 0.60966, 0.09233, 0.07102, 0.45147, 0.61214]
    b_0_0_c = [-0.37936, -0.19743, 0.22062, 0.21988, -0.15350, 0.01198]
    b_0_z_c = [0.30743, 0.21382, 0., 0., 0.27995, 0.21127]
    b_0_zz_c = [0., 0., -0.04419, -0.03749, 0., 0.]
    z_z_c = [0., 0., 0.78527, 0.92920, 0., 0.]
    b_V_0_c = [0.31475, 0.27075, -0.04805, -0.04629, 0.25471, 0.21535]
    b_V_z_c = [0.06073, 0.08202, -0.01454, -0.01833, 0.06761, 0.07763]
    V_SF_0 = dict(zip(cols, V_SF_0_c))
    V_SF_z = dict(zip(cols, V_SF_z_c))
    s_V_0 = dict(zip(cols, s_V_0_c))
    s_V_z = dict(zip(cols, s_V_z_c))
    b_0_0 = dict(zip(cols, b_0_0_c))
    b_0_z = dict(zip(cols, b_0_z_c))
    b_0_zz = dict(zip(cols, b_0_zz_c))
    z_z = dict(zip(cols, z_z_c))
    b_V_0 = dict(zip(cols, b_V_0_c))
    b_V_z = dict(zip(cols, b_V_z_c))

    # Compute the terms of the mass-dependant model at the given redshift
    #  All V's are in units of 220 km/s
    V_max_norm = V_max/220.
    V_SF = 10**(V_SF_0[case]+z*V_SF_z[case])
    s_V = s_V_0[case] + z*s_V_z[case]
    b_0 = b_0_0[case] + z*b_0_z[case]+b_0_zz[case]*(z-z_z[case])**2.
    b_V = b_V_0[case] + z*b_V_z[case]

    # Compute the terms of the final model
    b_x = 10**(0.5*(b_0+V_max_norm*b_V))
    s_o = s_V*abs(V_max_norm-V_SF)

    return([b_x, s_o, V_SF])


def b_s_Poole(s, V_max, z, halo_type, bias_type):

    """ This function expresses Equation (2) of Poole et al (2014)
        and fetches the parameters needed to compute it.

    Args:
        s (numpy.ndarray) : scale values

        V_max (float) :  halo maximum circular velocity

        z (float) :  redshift of interest

        halo_type (str) : halo type

        bias_type (str) : bias type

    Returns:
        A list containing two arrays with the values of `b_s`
        and `b_x` at each scale
    """

    # Set the bias parameters
    [b_x, s_o, V_SF] = set_b_s_Poole_params(V_max, z, halo_type, bias_type)

    # Create b(s) arrays
    V_max_norm = V_max/220.
    if(V_max_norm < V_SF):
        b_s = b_x*(1.-(s_o/s))**0.5
    else:
        b_s = b_x*(1.+(s_o/s))**0.5
    return([b_s, b_x])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
                                     'Bias model of Poole et al. (2014)')
    parser.add_argument('V_max', type=float, default=220.0,
                        help="maximum halo circular velocity")
    parser.add_argument('z', type=float, default=0.0,
                        help="redshift of interest")
    parser.add_argument('halo_type', type=str, default="FoF",
                        choices=["FoF", "sub"],
                        help="halo type")
    parser.add_argument('bias_type', type=str, default="zspace",
                        choices=["real", "zboost", "zspace"],
                        help="bias type")
    args = parser.parse_args()

    # Select the bias model you want
    V_max = args.V_max
    z = args.z
    halo_type = args.halo_type
    bias_type = args.bias_type

    # Print results
    print
    print "Selected mass      = %.1f km/s" % V_max
    print "Selected redshift  = %.2f" % z
    print "Selected halo_type = %s" % halo_type
    print "Selected bias_type = %s" % bias_type
    s = np.linspace(3, 150, 25)
    for i, s_i in enumerate(s):
        [b_s, b_x] = b_s_Poole(s_i, V_max, z, halo_type, bias_type)
        if(i == 0):
            print
            print "Large scale bias = %5.3f" % (b_x)
            print
            print "Scale dependant bias is:"
        print "%6.2f  %5.3f" % (s_i, b_s)

