#! /usr/bin/env python3
#

#%%
def normal_01_cdf_inv(cdf):

    # *****************************************************************************80
    #
    ## normal_01_cdf() evaluates the Normal 01 CDF.
    #
    #  Discussion:
    #
    #    The result is accurate to about 1 part in 10^16.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    28 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Michael Wichura,
    #    The Percentage Points of the Normal Distribution,
    #    Algorithm AS 241,
    #    Applied Statistics,
    #    Volume 37, Number 3, pages 477-484, 1988.
    #
    #  Input:
    #
    #    real CDF, the value of the CDF.
    #
    #  Output:
    #
    #    real VALUE, the argument of the CDF.
    #
    import numpy as np

    a = np.array(
        [
            3.3871328727963666080,
            1.3314166789178437745e2,
            1.9715909503065514427e3,
            1.3731693765509461125e4,
            4.5921953931549871457e4,
            6.7265770927008700853e4,
            3.3430575583588128105e4,
            2.5090809287301226727e3,
        ]
    )
    b = np.array(
        [
            1.0,
            4.2313330701600911252e1,
            6.8718700749205790830e2,
            5.3941960214247511077e3,
            2.1213794301586595867e4,
            3.9307895800092710610e4,
            2.8729085735721942674e4,
            5.2264952788528545610e3,
        ]
    )
    c = np.array(
        [
            1.42343711074968357734,
            4.63033784615654529590,
            5.76949722146069140550,
            3.64784832476320460504,
            1.27045825245236838258,
            2.41780725177450611770e-1,
            2.27238449892691845833e-2,
            7.74545014278341407640e-4,
        ]
    )
    const1 = 0.180625
    const2 = 1.6
    d = np.array(
        [
            1.0,
            2.05319162663775882187,
            1.67638483018380384940,
            6.89767334985100004550e-1,
            1.48103976427480074590e-1,
            1.51986665636164571966e-2,
            5.47593808499534494600e-4,
            1.05075007164441684324e-9,
        ]
    )
    e = np.array(
        [
            6.65790464350110377720,
            5.46378491116411436990,
            1.78482653991729133580,
            2.96560571828504891230e-1,
            2.65321895265761230930e-2,
            1.24266094738807843860e-3,
            2.71155556874348757815e-5,
            2.01033439929228813265e-7,
        ]
    )
    f = np.array(
        [
            1.0,
            5.99832206555887937690e-1,
            1.36929880922735805310e-1,
            1.48753612908506148525e-2,
            7.86869131145613259100e-4,
            1.84631831751005468180e-5,
            1.42151175831644588870e-7,
            2.04426310338993978564e-15,
        ]
    )
    huge = np.finfo(float).max
    split1 = 0.425
    split2 = 5.0

    if cdf <= 0.0:
        value = -huge
        return value

    if 1.0 <= cdf:
        value = huge
        return value

    q = cdf - 0.5

    if abs(q) <= split1:

        r = const1 - q * q
        value = q * r8poly_value_horner(7, a, r) / r8poly_value_horner(7, b, r)

    else:

        if q < 0.0:
            r = cdf
        else:
            r = 1.0 - cdf

        if r <= 0.0:

            value = huge

        else:

            r = np.sqrt(-np.log(r))

            if r <= split2:

                r = r - const2
                value = r8poly_value_horner(7, c, r) / r8poly_value_horner(7, d, r)

            else:

                r = r - split2
                value = r8poly_value_horner(7, e, r) / r8poly_value_horner(7, f, r)

        if q < 0.0:
            value = -value

    return value


def normal_01_cdf_inv_test():

    # *****************************************************************************80
    #
    ## normal_01_cdf_inv_test tests normal_01_cdf_inv.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    28 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_01_cdf_inv_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_cdf_inv inverts the CDF;")
    print("")
    print("      CDF             X                         X")
    print("                     (exact)                   (computed)")
    print("")

    n_data = 0

    while True:

        n_data, x1, cdf = normal_01_cdf_values(n_data)

        if n_data == 0:
            break

        x2 = normal_01_cdf_inv(cdf)

        print("  %14.6g  %24.16g  %24.16g" % (cdf, x1, x2))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_cdf_inv_test:")
    print("  Normal end of execution.")
    return


def normal_01_cdf(x):

    # *****************************************************************************80
    #
    ## normal_01_cdf evaluates the Normal 01 CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    27 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    A G Adams,
    #    Areas Under the Normal Curve,
    #    Algorithm 39,
    #    Computer j.,
    #    Volume 12, pages 197-198, 1969.
    #
    #  Input:
    #
    #    real X, the argument of the CDF.
    #
    #  Output:
    #
    #    real VALUE, the value of the CDF.
    #
    import numpy as np

    a1 = 0.398942280444
    a2 = 0.399903438504
    a3 = 5.75885480458
    a4 = 29.8213557808
    a5 = 2.62433121679
    a6 = 48.6959930692
    a7 = 5.92885724438
    b0 = 0.398942280385
    b1 = 3.8052e-08
    b2 = 1.00000615302
    b3 = 3.98064794e-04
    b4 = 1.98615381364
    b5 = 0.151679116635
    b6 = 5.29330324926
    b7 = 4.8385912808
    b8 = 15.1508972451
    b9 = 0.742380924027
    b10 = 30.789933034
    b11 = 3.99019417011
    #
    #  |X| <= 1.28.
    #
    if abs(x) <= 1.28:

        y = 0.5 * x * x

        q = 0.5 - abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))))
    #
    #  1.28 < |X| <= 12.7
    #
    elif abs(x) <= 12.7:

        y = 0.5 * x * x

        q = (
            np.exp(-y)
            * b0
            / (
                abs(x)
                - b1
                + b2
                / (
                    abs(x)
                    + b3
                    + b4
                    / (
                        abs(x)
                        - b5
                        + b6 / (abs(x) + b7 - b8 / (abs(x) + b9 + b10 / (abs(x) + b11)))
                    )
                )
            )
        )
    #
    #  12.7 < |X|
    #
    else:

        q = 0.0
    #
    #  Take account of negative X.
    #
    if x < 0.0:
        value = q
    else:
        value = 1.0 - q

    return value


def normal_01_cdf_test():

    # *****************************************************************************80
    #
    ## normal_01_cdf_test tests normal_01_cdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    27 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_01_cdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_cdf evaluates the CDF;")
    print("")
    print("       X              CDF                       CDF")
    print("                     (exact)                   (computed)")
    print("")

    n_data = 0

    while True:

        n_data, x, cdf1 = normal_01_cdf_values(n_data)

        if n_data == 0:
            break

        cdf2 = normal_01_cdf(x)

        print("  %14.6g  %24.16g  %24.16g" % (x, cdf1, cdf2))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_cdf_test:")
    print("  Normal end of execution.")
    return


def normal_01_cdf_values(n_data):

    # *****************************************************************************80
    #
    ## normal_01_cdf_values returns some values of the Normal 01 CDF.
    #
    #  Discussion:
    #
    #    In Mathematica, the function can be evaluated by:
    #
    #      Needs["Statistics`ContinuousDistributions`"]
    #      dist = NormalDistribution [ 0, 1 ]
    #      CDF [ dist, x ]
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    14 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Milton Abramowitz and Irene Stegun,
    #    Handbook of Mathematical Functions,
    #    US Department of Commerce, 1964.
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 17

    f_vec = np.array(
        (
            0.5000000000000000e00,
            0.5398278372770290e00,
            0.5792597094391030e00,
            0.6179114221889526e00,
            0.6554217416103242e00,
            0.6914624612740131e00,
            0.7257468822499270e00,
            0.7580363477769270e00,
            0.7881446014166033e00,
            0.8159398746532405e00,
            0.8413447460685429e00,
            0.9331927987311419e00,
            0.9772498680518208e00,
            0.9937903346742239e00,
            0.9986501019683699e00,
            0.9997673709209645e00,
            0.9999683287581669e00,
        )
    )

    x_vec = np.array(
        (
            0.0000000000000000e00,
            0.1000000000000000e00,
            0.2000000000000000e00,
            0.3000000000000000e00,
            0.4000000000000000e00,
            0.5000000000000000e00,
            0.6000000000000000e00,
            0.7000000000000000e00,
            0.8000000000000000e00,
            0.9000000000000000e00,
            0.1000000000000000e01,
            0.1500000000000000e01,
            0.2000000000000000e01,
            0.2500000000000000e01,
            0.3000000000000000e01,
            0.3500000000000000e01,
            0.4000000000000000e01,
        )
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        x = 0.0
        f = 0.0
    else:
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, x, f


def normal_01_cdf_values_test():

    # *****************************************************************************80
    #
    ## normal_01_cdf_values_test tests normal_01_cdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    14 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_01_cdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_cdf_values stores values of the unit normal CDF.")
    print("")
    print("      X         normal_01_cdf(X)")
    print("")

    n_data = 0

    while True:

        n_data, x, f = normal_01_cdf_values(n_data)

        if n_data == 0:
            break

        print("  %12f  %24.16f" % (x, f))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_cdf_values_test:")
    print("  Normal end of execution.")
    return


def normal_01_mean():

    # *****************************************************************************80
    #
    ## normal_01_mean returns the mean of the Normal 01 PDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Output:
    #
    #    real VALUE, the mean of the PDF.
    #
    value = 0.0

    return value


def normal_01_mean_test():

    # *****************************************************************************80
    #
    ## normal_01_mean_test tests normal_01_mean.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("normal_01_mean_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_mean computes the Normal 01 mean;")

    m = normal_01_mean()

    print("")
    print("  PDF mean =      %14g" % (m))

    nsample = 1000

    x = np.zeros(nsample)
    for i in range(0, nsample):
        x[i] = np.random.normal()

    print("")
    print("  Sample size =     %6d" % (nsample))
    print("  Sample mean =     %14g" % (np.mean(x)))
    print("  Sample maximum =  %14g" % (np.max(x)))
    print("  Sample minimum =  %14g" % (np.min(x)))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_mean_test:")
    print("  Normal end of execution.")
    return


def normal_01_moment(order):

    # *****************************************************************************80
    #
    ## normal_01_moment evaluates the moments of the Normal 01 PDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #
    #  Output:
    #
    #    real VALUE, the value of the moment.
    #
    from scipy.special import factorial2

    if (order % 2) == 0:
        value = factorial2(order - 1)
    else:
        value = 0.0

    return value


def normal_01_moment_test():

    # *****************************************************************************80
    #
    ## normal_01_moment_test tests normal_01_moment.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_01_moment_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_moment evaluates moments of the Normal 01 PDF;")
    print("")
    print("   Order     Moment")
    print("")

    for order in range(0, +11):

        moment = normal_01_moment(order)
        print("  %6d  %14.6g" % (order, moment))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_moment_test:")
    print("  Normal end of execution.")
    return


def normal_01_pdf(x):

    # *****************************************************************************80
    #
    ## normal_01_pdf evaluates the Normal 01 PDF.
    #
    #  Discussion:
    #
    #    The Normal 01 PDF is also called the "Standard Normal" PDF, or
    #    the Normal PDF with 0 mean and standard deviation 1.
    #
    #  Formula:
    #
    #    PDF(x) = exp ( - 0.5 * x^2 ) / sqrt ( 2 * pi )
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the PDF.
    #
    #  Output:
    #
    #    real VALUE, the value of the PDF.
    #
    import numpy as np

    value = np.exp(-0.5 * x * x) / np.sqrt(2.0 * np.pi)

    return value


def normal_01_pdf_test():

    # *****************************************************************************80
    #
    ## normal_01_pdf_test tests normal_01_pdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_01_pdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_pdf evaluates the PDF;")
    print("")
    print("       X              PDF")
    print("")

    for i in range(-20, +21):

        x = float(i) / 10.0
        pdf = normal_01_pdf(x)
        print("  %14.6g  %24.16g" % (x, pdf))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_pdf_test:")
    print("  Normal end of execution.")
    return


def normal_01_variance():

    # *****************************************************************************80
    #
    ## normal_01_variance returns the variance of the Normal 01 PDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Output:
    #
    #    real VALUE, the variance of the PDF.
    #
    value = 1.0

    return value


def normal_01_variance_test():

    # *****************************************************************************80
    #
    ## normal_01_variance_test tests normal_01_variance.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("normal_01_variance_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_01_variance computes the Normal 01 variance;")

    value = normal_01_variance()

    print("")
    print("  PDF variance =      %14g" % (value))

    nsample = 1000

    x = np.zeros(nsample)
    for i in range(0, nsample):
        x[i] = np.random.normal()

    value = r8vec_variance(nsample, x)

    print("")
    print("  Sample size =     %6d" % (nsample))
    print("  Sample variance = %14g" % (value))
    #
    #  Terminate.
    #
    print("")
    print("normal_01_variance_test:")
    print("  Normal end of execution.")
    return


def normal_ms_cdf_inv(cdf, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_cdf_inv inverts the CDF of the Normal MS distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real CDF, the value of the CDF.
    #    0.0 <= CDF <= 1.0.
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #    0.0 < SIGMA.
    #
    #  Output:
    #
    #    real VALUE, the corresponding argument.
    #
    if cdf < 0.0 or 1.0 < cdf:
        print("")
        print("normal_ms_cdf_inv - Fatal error!")
        print("  CDF < 0 or 1 < CDF.")
        raise Exception("normal_ms_cdf_inv - Fatal error!")

    y = normal_01_cdf_inv(cdf)

    value = mu + sigma * y

    return value


def normal_ms_cdf_inv_test():

    # *****************************************************************************80
    #
    ## normal_ms_cdf_inv_test tests normal_ms_cdf_inv.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_ms_cdf_inv_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_cdf_inv inverts the CDF")
    print("  of the Normal MS distribution.")

    mu = 100.0
    sigma = 15.0

    print("")
    print("  PDF parameter MU = %g" % (mu))
    print("  PDF parameter SIGMA = %g" % (sigma))

    print("")
    print("       X              CDF            CDF_inv")
    print("")

    for i in range(-20, +21):

        x = mu + sigma * float(i) / 10.0
        cdf = normal_ms_cdf(x, mu, sigma)
        x2 = normal_ms_cdf_inv(cdf, mu, sigma)
        print("  %14.6g  %14.6g  %14.6g" % (x, cdf, x2))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_cdf_inv_test:")
    print("  Normal end of execution.")
    return


def normal_ms_cdf(x, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_cdf() evaluates the CDF of the Normal MS distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the CDF.
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #    0.0 < SIGMA.
    #
    #  Output:
    #
    #    real VALUE, the value of the CDF.
    #
    y = (x - mu) / sigma

    value = normal_01_cdf(y)

    return value


def normal_ms_cdf_test():

    # *****************************************************************************80
    #
    ## normal_ms_cdf_test tests normal_ms_cdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_ms_cdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_cdf evaluates the CDF;")

    mu = 100.0
    sigma = 15.0

    print("")
    print("  PDF parameter MU = %g" % (mu))
    print("  PDF parameter SIGMA = %g" % (sigma))

    print("")
    print("       X              CDF")
    print("")

    for i in range(-20, +21):

        x = mu + sigma * float(i) / 10.0
        cdf = normal_ms_cdf(x, mu, sigma)
        print("  %14.6g  %24.16g" % (x, cdf))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_cdf_test:")
    print("  Normal end of execution.")
    return


def normal_ms_mean(mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_mean returns the mean of the Normal MS distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #    0.0 < SIGMA.
    #
    #  Output:
    #
    #    real VALUE, the mean of the PDF.
    #
    value = mu

    return value


def normal_ms_mean_test():

    # *****************************************************************************80
    #
    ## normal_ms_mean_test tests normal_ms_mean.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("normal_ms_mean_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_mean computes the mean")
    print("  of the Normal MS distribution.")

    mu = 100.0
    sigma = 15.0

    m = normal_ms_mean(mu, sigma)

    print("")
    print("  PDF parameter MU = %g" % (mu))
    print("  PDF parameter SIGMA = %g" % (sigma))
    print("  PDF mean = %g" % (m))

    nsample = 1000

    x = np.zeros(nsample)
    for i in range(0, nsample):
        x[i] = normal_ms_sample(mu, sigma)

    print("")
    print("  Sample size =     %6d" % (nsample))
    print("  Sample mean =     %14g" % (np.mean(x)))
    print("  Sample maximum =  %14g" % (np.max(x)))
    print("  Sample minimum =  %14g" % (np.min(x)))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_mean_test:")
    print("  Normal end of execution.")
    return


def normal_ms_moment_central(order, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_moment_central evaluates central moments of the Normal MS distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #  Output:
    #
    #    real VALUE, the value of the central moment.
    #
    from scipy.special import factorial2

    if (order % 2) == 0:
        value = factorial2(order - 1) * sigma**order
    else:
        value = 0.0

    return value


def normal_ms_moment_central_values(order, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_moment_central_values evaluates central moments 0 through 8 of the Normal PDF.
    #
    #  Discussion:
    #
    #    The formula was posted by John D Cook.
    #
    #    Order  Moment
    #    -----  ------
    #      0    1
    #      1    0
    #      2    sigma^2
    #      3    0
    #      4    3 sigma^4
    #      5    0
    #      6    15 sigma^6
    #      7    0
    #      8    105 sigma^8
    #      9    0
    #     10    945 sigma^10
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #    0 <= ORDER <= 8.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #  Output:
    #
    #    real VALUE, the value of the central moment.
    #
    if order == 0:
        value = 1.0
    elif order == 1:
        value = 0.0
    elif order == 2:
        value = sigma**2
    elif order == 3:
        value = 0.0
    elif order == 4:
        value = 3.0 * sigma**4
    elif order == 5:
        value = 0.0
    elif order == 6:
        value = 15.0 * sigma**6
    elif order == 7:
        value = 0.0
    elif order == 8:
        value = 105.0 * sigma**8
    elif order == 9:
        value = 0.0
    elif order == 10:
        value = 945.0 * sigma**10
    else:
        print("")
        print("normal_ms_moment_central_values - Fatal error!")
        print("  Only ORDERS 0 through 8 are available.")
        raise Exception("normal_ms_moment_central_values - Fatal error!")

    return value


def normal_ms_moment_central_test():

    # *****************************************************************************80
    #
    ## normal_ms_moment_central_test tests normal_ms_moment_central.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    test_num = 4
    mu_test = np.array([0.0, 2.0, 10.0, 0.0])
    sigma_test = np.array([1.0, 1.0, 2.0, 2.0])

    print("")
    print("normal_ms_moment_central_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_moment_central evaluates central moments")
    print("  of the Normal MS distribution.")

    for test in range(0, test_num):

        mu = mu_test[test]
        sigma = sigma_test[test]
        print("")
        print("  Mu = %g, Sigma = %g" % (mu, sigma))
        print(" Order  Moment")
        print("\n")

        for order in range(0, 9):
            moment1 = normal_ms_moment_central(order, mu, sigma)
            moment2 = normal_ms_moment_central_values(order, mu, sigma)
            print("  %2d  %12g  %12g" % (order, moment1, moment2))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_moment_central_test:")
    print("  Normal end of execution.")
    return


def normal_ms_moment(order, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_moment evaluates the moments of the Normal MS distribution.
    #
    #  Discussion:
    #
    #    The formula was posted by John D Cook.
    #
    #    Order  Moment
    #    -----  ------
    #      0    1
    #      1    mu
    #      2    mu ** 2 +         sigma ** 2
    #      3    mu ** 3 +  3 mu   sigma ** 2
    #      4    mu ** 4 +  6 mu ** 2 sigma ** 2 +   3      sigma ** 4
    #      5    mu ** 5 + 10 mu ** 3 sigma ** 2 +  15 mu   sigma ** 4
    #      6    mu ** 6 + 15 mu ** 4 sigma ** 2 +  45 mu ** 2 sigma ** 4 +  15      sigma ** 6
    #      7    mu ** 7 + 21 mu ** 5 sigma ** 2 + 105 mu ** 3 sigma ** 4 + 105 mu   sigma ** 6
    #      8    mu ** 8 + 28 mu ** 6 sigma ** 2 + 210 mu ** 4 sigma ** 4 + 420 mu ** 2 sigma ** 6 + 105 sigma ** 8
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #  Output:
    #
    #    real VALUE, the value of the moment.
    #
    from scipy.special import comb
    from scipy.special import factorial2

    j_hi = order // 2

    value = 0.0
    for j in range(0, j_hi + 1):
        value = value + comb(order, 2 * j) * factorial2(2 * j - 1) * mu ** (
            order - 2 * j
        ) * sigma ** (2 * j)

    return value


def normal_ms_moment_values(order, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_moment_values evaluates moments 0 through 8 of the Normal PDF.
    #
    #  Discussion:
    #
    #    The formula was posted by John D Cook.
    #
    #    Order  Moment
    #    -----  ------
    #      0    1
    #      1    mu
    #      2    mu ** 2 +         sigma ** 2
    #      3    mu ** 3 +  3 mu   sigma ** 2
    #      4    mu ** 4 +  6 mu ** 2 sigma ** 2 +   3      sigma ** 4
    #      5    mu ** 5 + 10 mu ** 3 sigma ** 2 +  15 mu   sigma ** 4
    #      6    mu ** 6 + 15 mu ** 4 sigma ** 2 +  45 mu ** 2 sigma ** 4 +  15      sigma ** 6
    #      7    mu ** 7 + 21 mu ** 5 sigma ** 2 + 105 mu ** 3 sigma ** 4 + 105 mu   sigma ** 6
    #      8    mu ** 8 + 28 mu ** 6 sigma ** 2 + 210 mu ** 4 sigma ** 4 + 420 mu ** 2 sigma ** 6 + 105 sigma ** 8
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #    0 <= ORDER <= 8.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #  Output:
    #
    #    real VALUE, the value of the central moment.
    #
    if order == 0:
        value = 1.0
    elif order == 1:
        value = mu
    elif order == 2:
        value = mu**2 + sigma**2
    elif order == 3:
        value = mu**3 + 3.0 * mu * sigma**2
    elif order == 4:
        value = mu**4 + 6.0 * mu**2 * sigma**2 + 3.0 * sigma**4
    elif order == 5:
        value = mu**5 + 10.0 * mu**3 * sigma**2 + 15.0 * mu * sigma**4
    elif order == 6:
        value = (
            mu**6
            + 15.0 * mu**4 * sigma**2
            + 45.0 * mu**2 * sigma**4
            + 15.0 * sigma**6
        )
    elif order == 7:
        value = (
            mu**7
            + 21.0 * mu**5 * sigma**2
            + 105.0 * mu**3 * sigma**4
            + 105.0 * mu * sigma**6
        )
    elif order == 8:
        value = (
            mu**8
            + 28.0 * mu**6 * sigma**2
            + 210.0 * mu**4 * sigma**4
            + 420.0 * mu**2 * sigma**6
            + 105.0 * sigma**8
        )
    else:
        print("")
        print("normal_ms_moment_values - Fatal error!")
        print("  Only ORDERS 0 through 8 are available.")
        raise Exception("normal_ms_moment_values - Fatal error!")

    return value


def normal_ms_moment_test():

    # *****************************************************************************80
    #
    ## normal_ms_moment_test tests normal_ms_moment.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    test_num = 4
    mu_test = np.array([0.0, 2.0, 10.0, 0.0])
    sigma_test = np.array([1.0, 1.0, 2.0, 2.0])

    print("")
    print("normal_ms_moment_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_moment evaluates moments of the Normal MS distribution.")

    for test in range(0, test_num):

        mu = mu_test[test]
        sigma = sigma_test[test]
        print("")
        print("  Mu = %g, Sigma = %g" % (mu, sigma))
        print(" Order  Moment")
        print("\n")

        for order in range(0, 9):
            moment1 = normal_ms_moment(order, mu, sigma)
            moment2 = normal_ms_moment_values(order, mu, sigma)
            print("  %2d  %12g  %12g" % (order, moment1, moment2))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_moment_test:")
    print("  Normal end of execution.")
    return


def normal_ms_pdf(x, mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_pdf evaluates the Normal MS PDF.
    #
    #  Discussion:
    #
    #    The Normal MS PDF is also called the Gaussian PDF.
    #
    #  Formula:
    #
    #    PDF(X)(MU,SIGMA) = EXP ( - 0.5 * ( ( X - MU ) / SIGMA )^2 )
    #      / SQRT ( 2 * PI * SIGMA^2 )
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the PDF.
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #    0.0 < SIGMA.
    #
    #  Output:
    #
    #    real VALUE, the value of the PDF.
    #
    import numpy as np

    value = np.exp(-0.5 * ((x - mu) / sigma) ** 2) / np.sqrt(2.0 * np.pi * sigma**2)

    return value


def normal_ms_pdf_test():

    # *****************************************************************************80
    #
    ## normal_ms_pdf_test tests normal_ms_pdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_ms_pdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_pdf evaluates the PDF")
    print("  for the Normal MS distribution.")

    mu = 100.0
    sigma = 15.0

    print("")
    print("  PDF parameter MU = %g" % (mu))
    print("  PDF parameter SIGMA = %g" % (sigma))

    print("")
    print("       X              PDF")
    print("")

    for i in range(-20, +21):

        x = mu + sigma * float(i) / 10.0
        pdf = normal_ms_pdf(x, mu, sigma)
        print("  %14.6g  %24.16g" % (x, pdf))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_pdf_test:")
    print("  Normal end of execution.")
    return


def normal_ms_sample(mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_sample samples the Normal MS distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #    0.0 < SIGMA.
    #
    #  Output:
    #
    #    real VALUE, a sample of the standard normal PDF.
    #
    import numpy as np

    y = np.random.normal()

    value = mu + sigma * y

    return value


def normal_ms_sample_test():

    # *****************************************************************************80
    #
    ## normal_ms_sample_test tests normal_ms_sample.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("normal_ms_sample_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_sample samples")
    print("  the Normal MS distribution.")

    mu = 100.0
    sigma = 15.0

    print("")
    print("  PDF parameter MU = %g\n" % (mu))
    print("  PDF parameter SIGMA = %g" % (sigma))

    print("")
    for i in range(0, 10):
        x = normal_ms_sample(mu, sigma)
        print("  %4d  %14.6g" % (i, x))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_sample_test")
    print("  Normal end of execution.")
    return


def normal_ms_variance(mu, sigma):

    # *****************************************************************************80
    #
    ## normal_ms_variance returns the variance of the Normal MS distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #    0.0 < SIGMA.
    #
    #  Output:
    #
    #    real VALUE, the variance of the PDF.
    #
    value = sigma * sigma

    return value


def normal_ms_variance_test():

    # *****************************************************************************80
    #
    ## normal_ms_variance_test tests normal_ms_variance.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("normal_ms_variance_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  normal_ms_variance computes the variance")
    print("  of the Normal MS distribution.")

    mu = 100.0
    sigma = 15.0

    value = normal_ms_variance(mu, sigma)

    print("")
    print("  PDF parameter MU = %g" % (mu))
    print("  PDF parameter SIGMA = %g" % (sigma))
    print("  PDF variance = %g" % (value))

    nsample = 1000

    x = np.zeros(nsample)
    for i in range(0, nsample):
        x[i] = normal_ms_sample(mu, sigma)

    value = r8vec_variance(nsample, x)

    print("")
    print("  Sample size =     %6d" % (nsample))
    print("  Sample variance = %14g" % (value))
    #
    #  Terminate.
    #
    print("")
    print("normal_ms_variance_test:")
    print("  Normal end of execution.")
    return


def r8_mop(i):

    # *****************************************************************************80
    #
    ## r8_mop returns the I-th power of -1 as an R8 value.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    01 June 2013
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer I, the power of -1.
    #
    #  Output:
    #
    #    real VALUE, the I-th power of -1.
    #
    if (i % 2) == 0:
        value = +1.0
    else:
        value = -1.0

    return value


def r8_mop_test():

    # *****************************************************************************80
    #
    ## r8_mop_test tests r8_mop.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    06 December 2014
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np

    print("")
    print("r8_mop_test")
    print("  r8_mop evaluates (-1.0)^I4 as an R8.")
    print("")
    print("    I4  r8_mop(I4)")
    print("")

    i4_min = -100
    i4_max = +100

    for test in range(0, 10):
        i4 = np.random.random_integers(i4_min, i4_max)
        r8 = r8_mop(i4)
        print("  %4d  %4.1f" % (i4, r8))
    #
    #  Terminate.
    #
    print("")
    print("r8_mop_test")
    print("  Normal end of execution.")
    return


def r8poly_print(m, a, title):

    # *****************************************************************************80
    #
    ## r8poly_print prints out a polynomial.
    #
    #  Discussion:
    #
    #    The power sum form is:
    #
    #      p(x) = a(0) + a(1) * x + ... + a(m-1) * x^(m-1) + a(m) * x^(m)
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    15 July 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer M, the nominal degree of the polynomial.
    #
    #    real A[0:M], the polynomial coefficients.
    #    A[0] is the constant term and
    #    A[M] is the coefficient of X^M.
    #
    #    string TITLE, a title.
    #
    if 0 < len(title):
        print("")
        print(title)
    print("")

    if a[m] < 0.0:
        plus_minus = "-"
    else:
        plus_minus = " "

    mag = abs(a[m])

    if 2 <= m:
        print("  p(x) = %c %g * x^%d" % (plus_minus, mag, m))
    elif m == 1:
        print("  p(x) = %c %g * x" % (plus_minus, mag))
    elif m == 0:
        print("  p(x) = %c %g" % (plus_minus, mag))

    for i in range(m - 1, -1, -1):

        if a[i] < 0.0:
            plus_minus = "-"
        else:
            plus_minus = "+"

        mag = abs(a[i])

        if mag != 0.0:

            if 2 <= i:
                print("         %c %g * x^%d" % (plus_minus, mag, i))
            elif i == 1:
                print("         %c %g * x" % (plus_minus, mag))
            elif i == 0:
                print("         %c %g" % (plus_minus, mag))


def r8poly_print_test():

    # *****************************************************************************80
    #
    ## r8poly_print_test tests r8poly_print.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 January 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("r8poly_print_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  r8poly_print prints an R8POLY.")

    m = 5
    c = np.array([12.0, -3.4, 56.0, 0.0, 0.78, 9.0])

    r8poly_print(m, c, "  The R8POLY:")
    #
    #  Terminate.
    #
    print("")
    print("r8poly_print_test:")
    print("  Normal end of execution.")

    return


def r8poly_value_horner(m, c, x):

    # *****************************************************************************80
    #
    ## r8poly_value_horner evaluates a polynomial using Horner's method.
    #
    #  Discussion:
    #
    #    The polynomial
    #
    #      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
    #
    #    is to be evaluated at the value X.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    05 January 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer M, the degree.
    #
    #    real C(0:M), the polynomial coefficients.
    #    C(I) is the coefficient of X^I.
    #
    #    real X, the evaluation point.
    #
    #  Output:
    #
    #    real VALUE, the polynomial value.
    #
    value = c[m]
    for i in range(m - 1, -1, -1):
        value = value * x + c[i]

    return value


def r8poly_value_horner_test():

    # *****************************************************************************80
    #
    ## r8poly_value_horner_test tests r8poly_value_horner.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    04 March 2016
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    m = 4
    n = 16
    c = np.array([24.0, -50.0, +35.0, -10.0, 1.0])

    print("")
    print("r8poly_value_horner_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  r8poly_value_horner evaluates a polynomial at a point")
    print("  using Horners method.")

    r8poly_print(m, c, "  The polynomial coefficients:")

    x_lo = 0.0
    x_hi = 5.0
    x = np.linspace(x_lo, x_hi, n)

    print("")
    print("   I    X    P(X)")
    print("")

    for i in range(0, n):
        p = r8poly_value_horner(m, c, x[i])
        print("  %2d  %8.4f  %14.6g" % (i, x[i], p))
    #
    #  Terminate.
    #
    print("")
    print("r8poly_value_horner_test:")
    print("  Normal end of execution.")
    return


def r8vec_linspace(n, a, b):

    # *****************************************************************************80
    #
    ## r8vec_linspace creates a column vector of linearly spaced values.
    #
    #  Discussion:
    #
    #    An R8VEC is a vector of R8's.
    #
    #    While MATLAB has the built in command
    #
    #      x = linspace ( a, b, n )
    #
    #    that command has the distinct disadvantage of returning a ROW vector.
    #
    #    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
    #
    #    In other words, the interval is divided into N-1 even subintervals,
    #    and the endpoints of intervals are used as the points.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    02 January 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer N, the number of entries in the vector.
    #
    #    real A, B, the first and last entries.
    #
    #  Output:
    #
    #    real X(N), a vector of linearly spaced data.
    #
    import numpy as np

    x = np.zeros(n)

    if n == 1:
        x[0] = (a + b) / 2.0
    else:
        for i in range(0, n):
            x[i] = ((n - 1 - i) * a + (i) * b) / (n - 1)

    return x


def r8vec_linspace_test():

    # *****************************************************************************80
    #
    ## r8vec_linspace_test tests r8vec_linspace.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    02 January 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("r8vec_linspace_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  r8vec_linspace returns evenly spaced values between A and B.")

    n = 5
    x_lo = 10.0
    x_hi = 20.0
    x = r8vec_linspace(n, x_lo, x_hi)

    r8vec_print(n, x, "  The linspace vector:")
    #
    #  Terminate.
    #
    print("")
    print("r8vec_linspace_test")
    print("  Normal end of execution.")
    return


def r8vec_print(n, a, title):

    # *****************************************************************************80
    #
    ## r8vec_print prints an R8VEC.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    31 August 2014
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer N, the dimension of the vector.
    #
    #    real A(N), the vector to be printed.
    #
    #    string TITLE, a title.
    #
    print("")
    print(title)
    print("")
    for i in range(0, n):
        print("%6d:  %12g" % (i, a[i]))


def r8vec_print_test():

    # *****************************************************************************80
    #
    ## r8vec_print_test tests r8vec_print.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    29 October 2014
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("r8vec_print_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  r8vec_print prints an R8VEC.")

    n = 4
    v = np.array([123.456, 0.000005, -1.0e06, 3.14159265], dtype=np.float64)
    r8vec_print(n, v, "  Here is an R8VEC:")
    #
    #  Terminate.
    #
    print("")
    print("r8vec_print_test:")
    print("  Normal end of execution.")
    return


def r8vec_variance(n, a):

    # *****************************************************************************80
    #
    ## r8vec_variance returns the variance of an R8VEC.
    #
    #  Discussion:
    #
    #    An R8VEC is a vector of R8's.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    03 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    integer N, the number of entries in the vector.
    #
    #    real A(N), the vector.
    #
    #  Output:
    #
    #    real VALUE, the variance of the vector.
    #
    import numpy as np

    #
    #  DDOF = 1 requests normalization by N-1 rather than N.
    #
    value = np.var(a, ddof=1)

    return value


def r8vec_variance_test():

    # *****************************************************************************80
    #
    ## r8vec_variance_test tests r8vec_variance.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    02 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    print("")
    print("r8vec_variance_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  r8vec_variance computes the variance of an R8VEC.")

    n = 10
    r8_lo = -5.0
    r8_hi = +5.0

    a = r8_lo + (r8_hi - r8_lo) * np.random.rand(n)

    r8vec_print(n, a, "  Input vector:")

    value = r8vec_variance(n, a)
    print("")
    print("  Value = %g" % (value))
    #
    #  Terminate.
    #
    print("")
    print("r8vec_variance_test:")
    print("  Normal end of execution.")
    return


def timestamp():

    # *****************************************************************************80
    #
    ## timestamp() prints the date as a timestamp.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    06 April 2013
    #
    #  Author:
    #
    #    John Burkardt
    #
    import time

    t = time.time()
    print(time.ctime(t))

    return None


def truncated_normal_ab_cdf_inv(cdf, mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf_inv inverts the truncated Normal CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real CDF, the value of the CDF.
    #    0.0 <= CDF <= 1.0.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #  Output:
    #
    #    real X, the corresponding argument.
    #
    if cdf < 0.0 or 1.0 < cdf:
        print("")
        print("truncated_normal_ab_cdf_inv - Fatal error!")
        print("  CDF < 0 or 1 < CDF.")
        raise Exception("truncated_normal_ab_cdf_inv - Fatal error!")

    alpha = (a - mu) / sigma
    beta = (b - mu) / sigma

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = normal_01_cdf(beta)

    xi_cdf = (beta_cdf - alpha_cdf) * cdf + alpha_cdf
    xi = normal_01_cdf_inv(xi_cdf)

    x = mu + sigma * xi

    return x


def truncated_normal_ab_cdf_inv_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf_inv_test tests truncated_normal_ab_cdf_inv.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    sample_num = 10
    a = 50.0
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_ab_cdf_inv_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_cdf_inv inverts the CDF of")
    print("  the Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,%g]" % (a, b))

    print("")
    print("             X            CDF            CDF_inv")
    print("")

    for i in range(0, sample_num):
        x = truncated_normal_ab_sample(mu, sigma, a, b)
        cdf = truncated_normal_ab_cdf(x, mu, sigma, a, b)
        x2 = truncated_normal_ab_cdf_inv(cdf, mu, sigma, a, b)
        print("  %14.6g  %14.6g  %14.6g" % (x, cdf, x2))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_cdf_inv_test")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_cdf(x, mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf evaluates the Truncated Normal CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    24 January 2017
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the CDF.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #  Output:
    #
    #    real VALUE, the value of the CDF.
    #
    if x < a:

        value = 0.0

    elif x <= b:

        alpha = (a - mu) / sigma
        beta = (b - mu) / sigma
        xi = (x - mu) / sigma

        alpha_cdf = normal_01_cdf(alpha)
        beta_cdf = normal_01_cdf(beta)
        xi_cdf = normal_01_cdf(xi)

        value = (xi_cdf - alpha_cdf) / (beta_cdf - alpha_cdf)

    else:

        value = 1.0

    return value


def truncated_normal_ab_cdf_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf_test tests truncated_normal_ab_cdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_ab_cdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_cdf evaluates the CDF")
    print("  of the Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean = mu")
    print("    standard deviation = sigma")
    print("  The parent distribution is truncated to")
    print("  the interval [a,b]")

    print("")
    print(
        "                                                           Stored         Computed"
    )
    print(
        "       X        Mu         S         A         B             CDF             CDF"
    )
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, b, x, cdf1 = truncated_normal_ab_cdf_values(n_data)

        if n_data == 0:
            break

        cdf2 = truncated_normal_ab_cdf(x, mu, sigma, a, b)

        print(
            "  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %14g  %14g"
            % (x, mu, sigma, a, b, cdf1, cdf2)
        )
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_cdf_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_cdf_values(n_data):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf_values: values of the Truncated Normal AB CDF.
    #
    #  Discussion:
    #
    #    The Normal distribution, with mean Mu and standard deviation Sigma,
    #    is truncated to the interval [A,B].
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 11

    a_vec = np.array((50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0))

    b_vec = np.array(
        (150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0)
    )

    f_vec = np.array(
        (
            0.3371694242213513,
            0.3685009225506048,
            0.4006444233448185,
            0.4334107066903040,
            0.4665988676496338,
            0.5000000000000000,
            0.5334011323503662,
            0.5665892933096960,
            0.5993555766551815,
            0.6314990774493952,
            0.6628305757786487,
        )
    )

    mu_vec = np.array(
        (100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0)
    )

    sigma_vec = np.array(
        (25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0)
    )

    x_vec = np.array(
        (90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 102.0, 104.0, 106.0, 108.0, 110.0)
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        mu = 0.0
        sigma = 0.0
        a = 0.0
        b = 0.0
        x = 0.0
        f = 0.0
    else:
        mu = mu_vec[n_data]
        sigma = sigma_vec[n_data]
        a = a_vec[n_data]
        b = b_vec[n_data]
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, mu, sigma, a, b, x, f


def truncated_normal_ab_cdf_values_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf_values_test tests truncated_normal_ab_cdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_ab_cdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print(
        "  truncated_normal_ab_cdf_values stores values of the truncated_normal_ab_cdf function."
    )
    print("")
    print(
        "            MU         SIGMA             A             B             X               F"
    )
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, b, x, f = truncated_normal_ab_cdf_values(n_data)

        if n_data == 0:
            break

        print("  %12g  %12g  %12g  %12g  %12g  %24.16g" % (mu, sigma, a, b, x, f))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_cdf_values_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_mean(mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_mean returns the mean of the Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the parent Normal Distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #  Output:
    #
    #    real VALUE, the mean of the PDF.
    #
    alpha = (a - mu) / sigma
    beta = (b - mu) / sigma

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = normal_01_cdf(beta)

    alpha_pdf = normal_01_pdf(alpha)
    beta_pdf = normal_01_pdf(beta)

    value = mu + sigma * (alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf)

    return value


def truncated_normal_ab_mean_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_mean_test tests truncated_normal_ab_mean.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    sample_num = 1000
    a = 50.0
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_ab_mean_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_mean computes the mean")
    print("  of the Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,%g]" % (a, b))

    m = truncated_normal_ab_mean(mu, sigma, a, b)

    print("")
    print("  PDF mean = %g" % (m))

    x = np.zeros(sample_num)
    for i in range(0, sample_num):
        x[i] = truncated_normal_ab_sample(mu, sigma, a, b)

    print("")
    print("  Sample size =     %6d" % (sample_num))
    print("  Sample mean =     %14g" % (np.mean(x)))
    print("  Sample maximum =  %14g" % (np.max(x)))
    print("  Sample minimum =  %14g" % (np.min(x)))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_mean_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_moment(order, mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_moment: moments of the truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Phoebus Dhrymes,
    #    Moments of Truncated Normal Distributions,
    #    May 2005.
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #    0 <= ORDER.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #    0 < S.
    #
    #    real A, B, the lower and upper truncation limits.
    #    A < B.
    #
    #  Output:
    #
    #    real VALUE, the moment of the PDF.
    #
    from scipy.special import comb

    if order < 0:
        print("")
        print("truncated_normal_ab_moment - Fatal error!")
        print("  ORDER < 0.")
        raise Exception("truncated_normal_ab_moment - Fatal error!")

    if sigma <= 0.0:
        print("")
        print("truncated_normal_ab_moment - Fatal error!")
        print("  SIGMA <= 0.0.")
        raise Exception("truncated_normal_ab_moment - Fatal error!")

    if b <= a:
        print("")
        print("truncated_normal_ab_moment - Fatal error!")
        print("  B <= A.")
        raise Exception("truncated_normal_ab_moment - Fatal error!")

    a_h = (a - mu) / sigma
    a_pdf = normal_01_pdf(a_h)
    a_cdf = normal_01_cdf(a_h)

    if a_cdf == 0.0:
        print("")
        print("truncated_normal_ab_moment - Fatal error!")
        print("  PDF/CDF ratio fails, because A_cdf is too small.")
        print("  A_pdf = %g" % (a_pdf))
        print("  A_cdf = %g" % (a_cdf))
        raise Exception("truncated_normal_ab_moment - Fatal error!")

    b_h = (b - mu) / sigma
    b_pdf = normal_01_pdf(b_h)
    b_cdf = normal_01_cdf(b_h)

    if b_cdf == 0.0:
        print("")
        print("truncated_normal_ab_moment - Fatal error!")
        print("  PDF/CDF ratio fails, because B_cdf too small.")
        print("  B_pdf = %g" % (b_pdf))
        print("  B_cdf = %g" % (b_cdf))
        raise Exception("truncated_normal_ab_moment - Fatal error!")

    value = 0.0
    irm2 = 0.0
    irm1 = 0.0

    for r in range(0, order + 1):

        if r == 0:
            ir = 1.0
        elif r == 1:
            ir = -(b_pdf - a_pdf) / (b_cdf - a_cdf)
        else:
            ir = (r - 1) * irm2 - (b_h ** (r - 1) * b_pdf - a_h ** (r - 1) * a_pdf) / (
                b_cdf - a_cdf
            )

        value = value + comb(order, r) * mu ** (order - r) * sigma**r * ir

        irm2 = irm1
        irm1 = ir

    return value


def truncated_normal_ab_moment_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_moment_test tests truncated_normal_ab_moment.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    test_num = 9
    mu_test = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 5.0])
    sigma_test = np.array([1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.5])
    a_test = np.array([-1.0, 0.0, -1.0, -1.0, 0.0, 0.5, -2.0, -4.0, 4.0])
    b_test = np.array([1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0, 4.0, 7.0])

    print("")
    print("truncated_normal_ab_moment_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_moment evaluates moments")
    print("  of the Truncated Normal distribution.")

    for test in range(0, test_num):

        mu = mu_test[test]
        sigma = sigma_test[test]
        a = a_test[test]
        b = b_test[test]
        print("")
        print(
            "  Test = %d, Mu = %g, Sigma = %g, A = %g, B = %g" % (test, mu, sigma, a, b)
        )
        print(" Order  Moment")
        print("\n")

        for order in range(0, 9):
            value = truncated_normal_ab_moment(order, mu, sigma, a, b)
            print("  %2d  %12g" % (order, value))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_moment_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_pdf(x, mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_pdf evaluates the Truncated Normal PDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    24 January 2017
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the PDF.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #  Output:
    #
    #    real VALUE, the value of the PDF.
    #
    if x < a:

        value = 0.0

    elif x <= b:

        alpha = (a - mu) / sigma
        beta = (b - mu) / sigma
        xi = (x - mu) / sigma

        alpha_cdf = normal_01_cdf(alpha)
        beta_cdf = normal_01_cdf(beta)
        xi_pdf = normal_01_pdf(xi)

        value = xi_pdf / (beta_cdf - alpha_cdf) / sigma

    else:

        value = 0.0

    return value


def truncated_normal_ab_pdf_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_pdf_test tests truncated_normal_ab_pdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_ab_pdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_pdf evaluates the PDF")
    print("  of the Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean = mu")
    print("    standard deviation = sigma")
    print("  The parent distribution is truncated to")
    print("  the interval [a,b]")

    print("")
    print(
        "                                                           Stored         Computed"
    )
    print(
        "       X        Mu         S         A         B             PDF             PDF"
    )
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, b, x, pdf1 = truncated_normal_ab_pdf_values(n_data)

        if n_data == 0:
            break

        pdf2 = truncated_normal_ab_pdf(x, mu, sigma, a, b)

        print(
            "  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %14g  %14g"
            % (x, mu, sigma, a, b, pdf1, pdf2)
        )
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_pdf_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_pdf_values(n_data):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_pdf_values: values of the Truncated Normal AB PDF.
    #
    #  Discussion:
    #
    #    The Normal distribution, with mean Mu and standard deviation Sigma,
    #    is truncated to the interval [A,B].
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 11

    a_vec = np.array((50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0))

    b_vec = np.array(
        (150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0)
    )

    f_vec = np.array(
        (
            0.01543301171801836,
            0.01588394472270638,
            0.01624375997031919,
            0.01650575046469259,
            0.01666496869385951,
            0.01671838200940538,
            0.01666496869385951,
            0.01650575046469259,
            0.01624375997031919,
            0.01588394472270638,
            0.01543301171801836,
        )
    )

    mu_vec = np.array(
        (100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0)
    )

    sigma_vec = np.array(
        (25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0)
    )

    x_vec = np.array(
        (90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 102.0, 104.0, 106.0, 108.0, 110.0)
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        mu = 0.0
        sigma = 0.0
        a = 0.0
        b = 0.0
        x = 0.0
        f = 0.0
    else:
        mu = mu_vec[n_data]
        sigma = sigma_vec[n_data]
        a = a_vec[n_data]
        b = b_vec[n_data]
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, mu, sigma, a, b, x, f


def truncated_normal_ab_pdf_values_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_pdf_values_test tests truncated_normal_ab_pdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_ab_pdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print(
        "  truncated_normal_ab_pdf_values stores values of the truncated_normal_ab_pdf function."
    )
    print("")
    print(
        "            MU         SIGMA             A             B             X               F"
    )
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, b, x, f = truncated_normal_ab_pdf_values(n_data)

        if n_data == 0:
            break

        print("  %12g  %12g  %12g  %12g  %12g  %24.16g" % (mu, sigma, a, b, x, f))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_pdf_values_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_sample(mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_sample samples the Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #  Output:
    #
    #    real X, a sample of the PDF.
    #
    import numpy as np

    alpha = (a - mu) / sigma
    beta = (b - mu) / sigma

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = normal_01_cdf(beta)

    u = np.random.rand()
    xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf)
    xi = normal_01_cdf_inv(xi_cdf)

    x = mu + sigma * xi

    return x


def truncated_normal_ab_sample_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_sample_test tests truncated_normal_ab_sample.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    sample_num = 10
    a = 50.0
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_ab_sample_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_sample samples")
    print("  the Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,%g]" % (a, b))
    print("")

    print("")
    for i in range(0, sample_num):
        x = truncated_normal_ab_sample(mu, sigma, a, b)
        print("  %4d  %14.6g" % (i, x))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_sample_test")
    print("  Normal end of execution.")
    return


def truncated_normal_ab_variance(mu, sigma, a, b):

    # *****************************************************************************80
    #
    ## truncated_normal_ab_variance: variance of the Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #
    #    real A, B, the lower and upper truncation limits.
    #
    #  Output:
    #
    #    real VALUE, the variance of the PDF.
    #
    alpha = (a - mu) / sigma
    beta = (b - mu) / sigma

    alpha_pdf = normal_01_pdf(alpha)
    beta_pdf = normal_01_pdf(beta)

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = normal_01_cdf(beta)

    value = (
        sigma
        * sigma
        * (
            1.0
            + (alpha * alpha_pdf - beta * beta_pdf) / (beta_cdf - alpha_cdf)
            - ((alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf)) ** 2
        )
    )

    return value


def truncated_normal_ab_variance_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_variance_test tests truncated_normal_ab_variance.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    sample_num = 1000
    a = 50.0
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_ab_variance_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_variance computes the variance")
    print("  of the Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,%g]" % (a, b))

    value = truncated_normal_ab_variance(mu, sigma, a, b)

    print("")
    print("  PDF variance = %g" % (value))

    x = np.zeros(sample_num)
    for i in range(0, sample_num):
        x[i] = truncated_normal_ab_sample(mu, sigma, a, b)

    value = r8vec_variance(sample_num, x)

    print("")
    print("  Sample size = %d" % (sample_num))
    print("  Sample variance = %g" % (value))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_variance_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_cdf_inv(cdf, mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_cdf_inv inverts the lower truncated Normal CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real CDF, the value of the CDF.
    #    0.0 <= CDF <= 1.0.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real X, the corresponding argument.
    #
    if cdf < 0.0 or 1.0 < cdf:
        print("")
        print("truncated_normal_a_cdf_inv - Fatal error!")
        print("  CDF < 0 or 1 < CDF.")
        raise Exception("truncated_normal_a_cdf_inv - Fatal error!")

    alpha = (a - mu) / sigma

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = 1.0

    xi_cdf = (beta_cdf - alpha_cdf) * cdf + alpha_cdf
    xi = normal_01_cdf_inv(xi_cdf)

    x = mu + sigma * xi

    return x


def truncated_normal_a_cdf_inv_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_cdf_inv_test tests truncated_normal_a_cdf_inv.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    sample_num = 10
    a = 50.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_a_cdf_inv_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_cdf_inv inverts the CDF of")
    print("  the lower Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,+oo)" % (a))

    print("")
    print("             X            CDF            CDF_inv")
    print("")

    for i in range(0, sample_num):
        x = truncated_normal_a_sample(mu, sigma, a)
        cdf = truncated_normal_a_cdf(x, mu, sigma, a)
        x2 = truncated_normal_a_cdf_inv(cdf, mu, sigma, a)
        print("  %14.6g  %14.6g  %14.6g" % (x, cdf, x2))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_cdf_inv_test")
    print("  Normal end of execution.")
    return


def truncated_normal_a_cdf(x, mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_cdf evaluates the lower Truncated Normal CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    24 January 2017
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the CDF.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the value of the CDF.
    #
    if x < a:

        value = 0.0

    else:

        alpha = (a - mu) / sigma
        xi = (x - mu) / sigma

        alpha_cdf = normal_01_cdf(alpha)
        beta_cdf = 1.0
        xi_cdf = normal_01_cdf(xi)

        value = (xi_cdf - alpha_cdf) / (beta_cdf - alpha_cdf)

    return value


def truncated_normal_a_cdf_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_cdf_test tests truncated_normal_a_cdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_a_cdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_cdf evaluates the CDF")
    print("  of the lower Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean = mu")
    print("    standard deviation = sigma")
    print("  The parent distribution is truncated to")
    print("  the interval [a,+oo)")

    print("")
    print("                                                 Stored         Computed")
    print("       X        Mu         S         A             CDF             CDF")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, x, cdf1 = truncated_normal_a_cdf_values(n_data)

        if n_data == 0:
            break

        cdf2 = truncated_normal_a_cdf(x, mu, sigma, a)

        print(
            "  %8.1f  %8.1f  %8.1f  %8.1f  %14g  %14g" % (x, mu, sigma, a, cdf1, cdf2)
        )
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_cdf_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_cdf_values(n_data):

    # *****************************************************************************80
    #
    ## truncated_normal_a_cdf_values: values of the Truncated Normal A CDF.
    #
    #  Discussion:
    #
    #    The Normal distribution, with mean Mu and standard deviation Sigma,
    #    is truncated to the interval [A,+oo).
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #    real A, the lower truncation limit.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 11

    a_vec = np.array((50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0))

    f_vec = np.array(
        (
            0.3293202045481688,
            0.3599223134505957,
            0.3913175216041539,
            0.4233210140873113,
            0.4557365629792204,
            0.4883601253415709,
            0.5209836877039214,
            0.5533992365958304,
            0.5854027290789878,
            0.6167979372325460,
            0.6474000461349729,
        )
    )

    mu_vec = np.array(
        (100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0)
    )

    sigma_vec = np.array(
        (25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0)
    )

    x_vec = np.array(
        (90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 102.0, 104.0, 106.0, 108.0, 110.0)
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        mu = 0.0
        sigma = 0.0
        a = 0.0
        x = 0.0
        f = 0.0
    else:
        mu = mu_vec[n_data]
        sigma = sigma_vec[n_data]
        a = a_vec[n_data]
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, mu, sigma, a, x, f


def truncated_normal_a_cdf_values_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_cdf_values_test tests truncated_normal_a_cdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_a_cdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print(
        "  truncated_normal_a_cdf_values stores values of the truncated_normal_a_cdf function."
    )
    print("")
    print("            MU         SIGMA             A             X               F")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, x, f = truncated_normal_a_cdf_values(n_data)

        if n_data == 0:
            break

        print("  %12g  %12g  %12g  %12g  %24.16g" % (mu, sigma, a, x, f))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_cdf_values_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_mean(mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_mean returns the mean of the lower Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the parent Normal Distribution.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the mean of the PDF.
    #
    alpha = (a - mu) / sigma

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = 1.0

    alpha_pdf = normal_01_pdf(alpha)
    beta_pdf = 0.0

    value = mu + sigma * (alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf)

    return value


def truncated_normal_a_mean_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_mean_test tests truncated_normal_a_mean.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    sample_num = 1000
    a = 50.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_a_mean_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_mean computes the mean")
    print("  of the Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,+oo)" % (a))

    m = truncated_normal_a_mean(mu, sigma, a)

    print("")
    print("  PDF mean = %g" % (m))

    x = np.zeros(sample_num)
    for i in range(0, sample_num):
        x[i] = truncated_normal_a_sample(mu, sigma, a)

    print("")
    print("  Sample size =     %6d" % (sample_num))
    print("  Sample mean =     %14g" % (np.mean(x)))
    print("  Sample maximum =  %14g" % (np.max(x)))
    print("  Sample minimum =  %14g" % (np.min(x)))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_mean_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_moment(order, mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_moment: moments of the truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Phoebus Dhrymes,
    #    Moments of Truncated Normal Distributions,
    #    May 2005.
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #    0 <= ORDER.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #    0 < S.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the moment of the PDF.
    #
    value = r8_mop(order) * truncated_normal_b_moment(order, -mu, sigma, -a)

    return value


def truncated_normal_a_moment_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_moment_test tests truncated_normal_a_moment.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    test_num = 6
    mu_test = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -5.0])
    sigma_test = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 1.0])
    a_test = np.array([0.0, -10.0, 10.0, -10.0, 10.0, -10.0])

    print("")
    print("truncated_normal_a_moment_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_moment evaluates moments")
    print("  of the lower Truncated Normal distribution.")

    for test in range(0, test_num):

        mu = mu_test[test]
        sigma = sigma_test[test]
        a = a_test[test]
        print("")
        print("  Test = %d, Mu = %g, Sigma = %g, A = %g" % (test, mu, sigma, a))
        print(" Order  Moment")
        print("\n")

        for order in range(0, 9):
            value = truncated_normal_a_moment(order, mu, sigma, a)
            print("  %2d  %12g" % (order, value))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_moment_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_pdf(x, mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_pdf evaluates the lower Truncated Normal PDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    24 January 2017
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the PDF.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the value of the PDF.
    #
    if x < a:

        value = 0.0

    else:

        alpha = (a - mu) / sigma
        xi = (x - mu) / sigma

        alpha_cdf = normal_01_cdf(alpha)
        beta_cdf = 1.0
        xi_pdf = normal_01_pdf(xi)

        value = xi_pdf / (beta_cdf - alpha_cdf) / sigma

    return value


def truncated_normal_a_pdf_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_pdf_test tests truncated_normal_a_pdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_a_pdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_pdf evaluates the PDF")
    print("  of the lower Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean = mu")
    print("    standard deviation = sigma")
    print("  The parent distribution is truncated to")
    print("  the interval [a,+oo)")

    print("")
    print("                                                 Stored         Computed")
    print("       X        Mu         S         A             PDF             PDF")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, x, pdf1 = truncated_normal_a_pdf_values(n_data)

        if n_data == 0:
            break

        pdf2 = truncated_normal_a_pdf(x, mu, sigma, a)

        print(
            "  %8.1f  %8.1f  %8.1f  %8.1f  %14g  %14g" % (x, mu, sigma, a, pdf1, pdf2)
        )
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_pdf_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_pdf_values(n_data):

    # *****************************************************************************80
    #
    ## truncated_normal_a_pdf_values: values of the Truncated Normal A PDF.
    #
    #  Discussion:
    #
    #    The Normal distribution, with mean Mu and standard deviation Sigma,
    #    is truncated to the interval [A,+oo).
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #    real A, the lower truncation limit.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 11

    a_vec = np.array((50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0))

    f_vec = np.array(
        (
            0.01507373507401876,
            0.01551417047139894,
            0.01586560931024694,
            0.01612150073158793,
            0.01627701240029317,
            0.01632918226724295,
            0.01627701240029317,
            0.01612150073158793,
            0.01586560931024694,
            0.01551417047139894,
            0.01507373507401876,
        )
    )

    mu_vec = np.array(
        (100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0)
    )

    sigma_vec = np.array(
        (25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0)
    )

    x_vec = np.array(
        (90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 102.0, 104.0, 106.0, 108.0, 110.0)
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        mu = 0.0
        sigma = 0.0
        a = 0.0
        x = 0.0
        f = 0.0
    else:
        mu = mu_vec[n_data]
        sigma = sigma_vec[n_data]
        a = a_vec[n_data]
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, mu, sigma, a, x, f


def truncated_normal_a_pdf_values_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_pdf_values_test tests truncated_normal_a_pdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_a_pdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print(
        "  truncated_normal_a_pdf_values stores values of the truncated_normal_a_pdf function."
    )
    print("")
    print("            MU         SIGMA             A             X               F")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, a, x, f = truncated_normal_a_pdf_values(n_data)

        if n_data == 0:
            break

        print("  %12g  %12g  %12g  %12g  %24.16g" % (mu, sigma, a, x, f))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_pdf_values_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_a_sample(mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_sample samples the lower Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real X, a sample of the PDF.
    #
    import numpy as np

    alpha = (a - mu) / sigma

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = 1.0

    u = np.random.rand()
    xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf)
    xi = normal_01_cdf_inv(xi_cdf)

    x = mu + sigma * xi

    return x


def truncated_normal_a_sample_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_sample_test tests truncated_normal_a_sample.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    sample_num = 10
    a = 50.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_a_sample_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_sample samples")
    print("  the lower Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,+oo)" % (a))
    print("")

    print("")
    for i in range(0, sample_num):
        x = truncated_normal_a_sample(mu, sigma, a)
        print("  %4d  %14.6g" % (i, x))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_sample_test")
    print("  Normal end of execution.")
    return


def truncated_normal_a_variance(mu, sigma, a):

    # *****************************************************************************80
    #
    ## truncated_normal_a_variance: variance of the lower Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #
    #    real A, the lower truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the variance of the PDF.
    #
    alpha = (a - mu) / sigma

    alpha_pdf = normal_01_pdf(alpha)
    beta_pdf = 0.0

    alpha_cdf = normal_01_cdf(alpha)
    beta_cdf = 1.0

    value = (
        sigma
        * sigma
        * (
            1.0
            + (alpha * alpha_pdf - 0.0) / (beta_cdf - alpha_cdf)
            - ((alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf)) ** 2
        )
    )

    return value


def truncated_normal_a_variance_test():

    # *****************************************************************************80
    #
    ## truncated_normal_a_variance_test tests truncated_normal_a_variance.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    sample_num = 1000
    a = 50.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_a_variance_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_a_variance computes the variance")
    print("  of the Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval [%g,+oo)" % (a))

    value = truncated_normal_a_variance(mu, sigma, a)

    print("")
    print("  PDF variance = %g" % (value))

    x = np.zeros(sample_num)
    for i in range(0, sample_num):
        x[i] = truncated_normal_a_sample(mu, sigma, a)

    value = r8vec_variance(sample_num, x)

    print("")
    print("  Sample size = %d" % (sample_num))
    print("  Sample variance = %g" % (value))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_a_variance_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_cdf_inv(cdf, mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_cdf_inv inverts the upper truncated Normal CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real CDF, the value of the CDF.
    #    0.0 <= CDF <= 1.0.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real B, the upper truncation limit.
    #
    #  Output:
    #
    #    real X, the corresponding argument.
    #
    if cdf < 0.0 or 1.0 < cdf:
        print("")
        print("truncated_normal_b_cdf_inv - Fatal error!")
        print("  CDF < 0 or 1 < CDF.")
        raise Exception("truncated_normal_b_cdf_inv - Fatal error!")

    beta = (b - mu) / sigma

    alpha_cdf = 0.0
    beta_cdf = normal_01_cdf(beta)

    xi_cdf = (beta_cdf - alpha_cdf) * cdf + alpha_cdf
    xi = normal_01_cdf_inv(xi_cdf)

    x = mu + sigma * xi

    return x


def truncated_normal_b_cdf_inv_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_cdf_inv_test tests truncated_normal_b_cdf_inv.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    sample_num = 10
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_b_cdf_inv_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_b_cdf_inv inverts the CDF of")
    print("  the Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval (-oo,%g]" % (b))

    print("")
    print("             X            CDF            CDF_inv")
    print("")

    for i in range(0, sample_num):
        x = truncated_normal_b_sample(mu, sigma, b)
        cdf = truncated_normal_b_cdf(x, mu, sigma, b)
        x2 = truncated_normal_b_cdf_inv(cdf, mu, sigma, b)
        print("  %14.6g  %14.6g  %14.6g" % (x, cdf, x2))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_cdf_inv_test")
    print("  Normal end of execution.")
    return


def truncated_normal_b_cdf(x, mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_cdf evaluates the upper Truncated Normal CDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    24 January 2017
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the CDF.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real B, the upper truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the value of the CDF.
    #
    if x <= b:

        beta = (b - mu) / sigma
        xi = (x - mu) / sigma

        alpha_cdf = 0.0
        beta_cdf = normal_01_cdf(beta)
        xi_cdf = normal_01_cdf(xi)

        value = (xi_cdf - alpha_cdf) / (beta_cdf - alpha_cdf)

    else:

        value = 1.0

    return value


def truncated_normal_b_cdf_test():

    # *****************************************************************************80
    #
    ## truncated_normal_ab_cdf_test tests truncated_normal_ab_cdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_b_cdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_ab_cdf evaluates the CDF")
    print("  of the upper Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean = mu")
    print("    standard deviation = sigma")
    print("  The parent distribution is truncated to")
    print("  the interval (-oo,b]")

    print("")
    print("                                                 Stored         Computed")
    print("       X        Mu         S         B             CDF             CDF")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, b, x, cdf1 = truncated_normal_b_cdf_values(n_data)

        if n_data == 0:
            break

        cdf2 = truncated_normal_b_cdf(x, mu, sigma, b)

        print(
            "  %8.1f  %8.1f  %8.1f  %8.1f  %14g  %14g" % (x, mu, sigma, b, cdf1, cdf2)
        )
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_cdf_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_cdf_values(n_data):

    # *****************************************************************************80
    #
    ## truncated_normal_b_cdf_values: values of the Truncated Normal B CDF.
    #
    #  Discussion:
    #
    #    The Normal distribution, with mean Mu and standard deviation Sigma,
    #    is truncated to the interval (-oo,B].
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #    real B, the upper truncation limit.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 11

    b_vec = np.array(
        (150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0)
    )

    f_vec = np.array(
        (
            0.3525999538650271,
            0.3832020627674540,
            0.4145972709210122,
            0.4466007634041696,
            0.4790163122960786,
            0.5116398746584291,
            0.5442634370207796,
            0.5766789859126887,
            0.6086824783958461,
            0.6400776865494043,
            0.6706797954518312,
        )
    )

    mu_vec = np.array(
        (100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0)
    )

    sigma_vec = np.array(
        (25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0)
    )

    x_vec = np.array(
        (90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 102.0, 104.0, 106.0, 108.0, 110.0)
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        mu = 0.0
        sigma = 0.0
        b = 0.0
        x = 0.0
        f = 0.0
    else:
        mu = mu_vec[n_data]
        sigma = sigma_vec[n_data]
        b = b_vec[n_data]
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, mu, sigma, b, x, f


def truncated_normal_b_cdf_values_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_cdf_values_test tests truncated_normal_b_cdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_b_cdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print(
        "  truncated_normal_b_cdf_values stores values of the truncated_normal_b_cdf function."
    )
    print("")
    print("            MU         SIGMA             B             X               F")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, b, x, f = truncated_normal_b_cdf_values(n_data)

        if n_data == 0:
            break

        print("  %12g  %12g  %12g  %12g  %24.16g" % (mu, sigma, b, x, f))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_cdf_values_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_mean(mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_mean returns the mean of the upper Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the parent Normal Distribution.
    #
    #    real B, the upper truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the mean of the PDF.
    #
    beta = (b - mu) / sigma

    alpha_cdf = 0.0
    beta_cdf = normal_01_cdf(beta)

    alpha_pdf = 0.0
    beta_pdf = normal_01_pdf(beta)

    value = mu + sigma * (alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf)

    return value


def truncated_normal_b_mean_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_mean_test tests truncated_normal_b_mean.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    sample_num = 1000
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_b_mean_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_b_mean computes the mean")
    print("  of the Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval (-oo,%g]" % (b))

    m = truncated_normal_b_mean(mu, sigma, b)

    print("")
    print("  PDF mean = %g" % (m))

    x = np.zeros(sample_num)
    for i in range(0, sample_num):
        x[i] = truncated_normal_b_sample(mu, sigma, b)

    print("")
    print("  Sample size =     %6d" % (sample_num))
    print("  Sample mean =     %14g" % (np.mean(x)))
    print("  Sample maximum =  %14g" % (np.max(x)))
    print("  Sample minimum =  %14g" % (np.min(x)))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_mean_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_moment(order, mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_moment: moments of upper truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Phoebus Dhrymes,
    #    Moments of Truncated Normal Distributions,
    #    May 2005.
    #
    #  Input:
    #
    #    integer ORDER, the order of the moment.
    #    0 <= ORDER.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #    0 < S.
    #
    #    real B, the upper truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the moment of the PDF.
    #
    from scipy.special import comb

    if order < 0:
        print("")
        print("truncated_normal_b_moment - Fatal error!")
        print("  ORDER < 0.")
        raise Exception("truncated_normal_b_moment - Fatal error!")

    if sigma <= 0.0:
        print("")
        print("truncated_normal_b_moment - Fatal error!")
        print("  SIGMA <= 0.0.")
        raise Exception("truncated_normal_b_moment - Fatal error!")

    b_h = (b - mu) / sigma
    b_pdf = normal_01_pdf(b_h)
    b_cdf = normal_01_cdf(b_h)

    if b_cdf == 0.0:
        print("")
        print("truncated_normal_b_moment - Fatal error!")
        print("  PDF/CDF ratio fails, because B_cdf too small.")
        print("  B_pdf = %g" % (b_pdf))
        print("  B_cdf = %g" % (b_cdf))
        raise Exception("truncated_normal_b_moment - Fatal error!")

    f = b_pdf / b_cdf

    value = 0.0
    irm2 = 0.0
    irm1 = 0.0

    for r in range(0, order + 1):

        if r == 0:
            ir = 1.0
        elif r == 1:
            ir = -f
        else:
            ir = -(b_h ** (r - 1)) * f + (r - 1) * irm2

        value = value + comb(order, r) * mu ** (order - r) * sigma**r * ir

        irm2 = irm1
        irm1 = ir

    return value


def truncated_normal_b_moment_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_moment_test tests truncated_normal_b_moment.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    test_num = 6
    mu_test = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 5.0])
    sigma_test = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 1.0])
    b_test = np.array([0.0, 10.0, -10.0, 10.0, -10.0, 10.0])

    print("")
    print("truncated_normal_b_moment_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_b_moment evaluates moments")
    print("  of the upper Truncated Normal distribution.")

    for test in range(0, test_num):

        mu = mu_test[test]
        sigma = sigma_test[test]
        b = b_test[test]
        print("")
        print("  Test = %d, Mu = %g, Sigma = %g, B = %g" % (test, mu, sigma, b))
        print(" Order  Moment")
        print("\n")

        for order in range(0, 9):
            value = truncated_normal_b_moment(order, mu, sigma, b)
            print("  %2d  %12g" % (order, value))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_moment_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_pdf(x, mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_pdf evaluates the upper Truncated Normal PDF.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    24 January 2017
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real X, the argument of the PDF.
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real B, the upper truncation limit.
    #
    #  Output:
    #
    #    real VALUE, the value of the PDF.
    #
    if x <= b:

        beta = (b - mu) / sigma
        xi = (x - mu) / sigma

        alpha_cdf = 0.0
        beta_cdf = normal_01_cdf(beta)
        xi_pdf = normal_01_pdf(xi)

        value = xi_pdf / (beta_cdf - alpha_cdf) / sigma

    else:

        value = 0.0

    return value


def truncated_normal_b_pdf_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_pdf_test tests truncated_normal_b_pdf.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_b_pdf_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_b_pdf evaluates the PDF")
    print("  of the upper Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean = mu")
    print("    standard deviation = sigma")
    print("  The parent distribution is truncated to")
    print("  the interval (-oo,b]")

    print("")
    print("                                                 Stored         Computed")
    print("       X        Mu         S         B             PDF             PDF")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, b, x, pdf1 = truncated_normal_b_pdf_values(n_data)

        if n_data == 0:
            break

        pdf2 = truncated_normal_b_pdf(x, mu, sigma, b)

        print(
            "  %8.1f  %8.1f  %8.1f  %8.1f  %14g  %14g" % (x, mu, sigma, b, pdf1, pdf2)
        )
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_ab_pdf_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_pdf_values(n_data):

    # *****************************************************************************80
    #
    ## truncated_normal_b_pdf_values: values of the Truncated Normal B PDF.
    #
    #  Discussion:
    #
    #    The Normal distribution, with mean Mu and standard deviation Sigma,
    #    is truncated to the interval (-oo,B].
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Stephen Wolfram,
    #    The Mathematica Book,
    #    Fourth Edition,
    #    Wolfram Media / Cambridge University Press, 1999.
    #
    #  Input:
    #
    #    integer N_DATA.  The user sets N_DATA to 0 before the first call.
    #
    #  Output:
    #
    #    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
    #    returns the corresponding data; when there is no more data, the
    #    output value of N_DATA will be 0 again.
    #
    #    real MU, the mean of the distribution.
    #
    #    real SIGMA, the standard deviation of the distribution.
    #
    #    real B, the upper truncation limit.
    #
    #    real X, the argument of the function.
    #
    #    real F, the value of the function.
    #
    import numpy as np

    n_max = 11

    b_vec = np.array(
        (150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0)
    )

    f_vec = np.array(
        (
            0.01507373507401876,
            0.01551417047139894,
            0.01586560931024694,
            0.01612150073158793,
            0.01627701240029317,
            0.01632918226724295,
            0.01627701240029317,
            0.01612150073158793,
            0.01586560931024694,
            0.01551417047139894,
            0.01507373507401876,
        )
    )

    mu_vec = np.array(
        (100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0)
    )

    sigma_vec = np.array(
        (25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0)
    )

    x_vec = np.array(
        (90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 102.0, 104.0, 106.0, 108.0, 110.0)
    )

    if n_data < 0:
        n_data = 0

    if n_max <= n_data:
        n_data = 0
        mu = 0.0
        sigma = 0.0
        b = 0.0
        x = 0.0
        f = 0.0
    else:
        mu = mu_vec[n_data]
        sigma = sigma_vec[n_data]
        b = b_vec[n_data]
        x = x_vec[n_data]
        f = f_vec[n_data]
        n_data = n_data + 1

    return n_data, mu, sigma, b, x, f


def truncated_normal_b_pdf_values_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_pdf_values_test tests truncated_normal_b_pdf_values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    22 February 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_b_pdf_values_test:")
    print("  Python version: %s" % (platform.python_version()))
    print(
        "  truncated_normal_b_pdf_values stores values of the truncated_normal_b_pdf function."
    )
    print("")
    print("            MU         SIGMA             B             X               F")
    print("")

    n_data = 0

    while True:

        n_data, mu, sigma, b, x, f = truncated_normal_b_pdf_values(n_data)

        if n_data == 0:
            break

        print("  %12g  %12g  %12g  %12g  %24.16g" % (mu, sigma, b, x, f))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_pdf_values_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_b_sample(mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_sample samples the upper Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the mean and standard deviation of the
    #    parent Normal distribution.
    #
    #    real B, the upper truncation limit.
    #
    #  Output:
    #
    #    real X, a sample of the PDF.
    #
    import numpy as np

    beta = (b - mu) / sigma

    alpha_cdf = 0.0
    beta_cdf = normal_01_cdf(beta)

    u = np.random.rand()
    xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf)
    xi = normal_01_cdf_inv(xi_cdf)

    x = mu + sigma * xi

    return x


def truncated_normal_b_sample_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_sample_test tests truncated_normal_b_sample.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    sample_num = 10
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_b_sample_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_b_sample samples")
    print("  the upper Truncated Normal distribution.")

    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval (-oo,%g]" % (b))
    print("")

    print("")
    for i in range(0, sample_num):
        x = truncated_normal_b_sample(mu, sigma, b)
        print("  %4d  %14.6g" % (i, x))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_sample_test")
    print("  Normal end of execution.")
    return


def truncated_normal_b_variance(mu, sigma, b):

    # *****************************************************************************80
    #
    ## truncated_normal_b_variance: variance of the upper Truncated Normal distribution.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Input:
    #
    #    real MU, SIGMA, the parameters of the PDF.
    #
    #    real B, the upper truncation limits.
    #
    #  Output:
    #
    #    real VALUE, the variance of the PDF.
    #
    beta = (b - mu) / sigma

    alpha_pdf = 0.0
    beta_pdf = normal_01_pdf(beta)

    alpha_cdf = 0.0
    beta_cdf = normal_01_cdf(beta)

    value = (
        sigma
        * sigma
        * (
            1.0
            + (0.0 - beta * beta_pdf) / (beta_cdf - alpha_cdf)
            - ((alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf)) ** 2
        )
    )

    return value


def truncated_normal_b_variance_test():

    # *****************************************************************************80
    #
    ## truncated_normal_b_variance_test tests truncated_normal_b_variance.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import numpy as np
    import platform

    sample_num = 1000
    b = 150.0
    mu = 100.0
    sigma = 25.0

    print("")
    print("truncated_normal_b_variance_test")
    print("  Python version: %s" % (platform.python_version()))
    print("  truncated_normal_b_variance computes the variance")
    print("  of the Truncated Normal distribution.")
    print("")
    print('  The "parent" normal distribution has')
    print("    mean =               %g" % (mu))
    print("    standard deviation = %g" % (sigma))
    print("  The parent distribution is truncated to")
    print("  the interval (-oo,%g]" % (b))

    value = truncated_normal_b_variance(mu, sigma, b)

    print("")
    print("  PDF variance = %g" % (value))

    x = np.zeros(sample_num)
    for i in range(0, sample_num):
        x[i] = truncated_normal_b_sample(mu, sigma, b)

    value = r8vec_variance(sample_num, x)

    print("")
    print("  Sample size = %d" % (sample_num))
    print("  Sample variance = %g" % (value))
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_b_variance_test:")
    print("  Normal end of execution.")
    return


def truncated_normal_test():

    # *****************************************************************************80
    #
    ## truncated_normal_test() tests truncated_normal().
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    09 March 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    import platform

    print("")
    print("truncated_normal_test()")
    print("  Python version: %s" % (platform.python_version()))
    print("  Test truncated_normal().")
    #
    #  Utilities:
    #
    normal_01_cdf_values_test()

    r8_mop_test()

    r8poly_print_test()
    r8poly_value_horner_test()

    r8vec_linspace_test()
    r8vec_print_test()
    r8vec_variance_test()

    truncated_normal_a_cdf_values_test()
    truncated_normal_a_pdf_values_test()

    truncated_normal_ab_cdf_values_test()
    truncated_normal_ab_pdf_values_test()

    truncated_normal_b_cdf_values_test()
    truncated_normal_b_pdf_values_test()
    #
    #  Library:
    #
    normal_01_cdf_test()
    normal_01_cdf_inv_test()
    normal_01_mean_test()
    normal_01_moment_test()
    normal_01_pdf_test()
    normal_01_variance_test()

    normal_ms_cdf_test()
    normal_ms_cdf_inv_test()
    normal_ms_mean_test()
    normal_ms_moment_test()
    normal_ms_moment_central_test()
    normal_ms_pdf_test()
    normal_ms_sample_test()
    normal_ms_variance_test()

    truncated_normal_a_cdf_test()
    truncated_normal_a_cdf_inv_test()
    truncated_normal_a_mean_test()
    truncated_normal_a_moment_test()
    truncated_normal_a_pdf_test()
    truncated_normal_a_sample_test()
    truncated_normal_a_variance_test()

    truncated_normal_ab_cdf_test()
    truncated_normal_ab_cdf_inv_test()
    truncated_normal_ab_mean_test()
    truncated_normal_ab_moment_test()
    truncated_normal_ab_pdf_test()
    truncated_normal_ab_sample_test()
    truncated_normal_ab_variance_test()

    truncated_normal_b_cdf_test()
    truncated_normal_b_cdf_inv_test()
    truncated_normal_b_mean_test()
    truncated_normal_b_moment_test()
    truncated_normal_b_pdf_test()
    truncated_normal_b_sample_test()
    truncated_normal_b_variance_test()
    #
    #  Terminate.
    #
    print("")
    print("truncated_normal_test():")
    print("  Normal end of execution.")
    return


if __name__ == "__main__":
    timestamp()
    truncated_normal_test()
    timestamp()

# %%
