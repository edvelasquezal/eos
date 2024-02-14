/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2024 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2018 Frederik Beaujean
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_MATHS_INTEGRATE_IMPL_HH
#define EOS_GUARD_EOS_MATHS_INTEGRATE_IMPL_HH 1

#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-cubature.hh>
#include <eos/maths/matrix.hh>
#include <eos/utils/log.hh>

#include <cassert>
#include <numeric>
#include <vector>

#include <iostream>

namespace eos
{
    namespace custom
    {
        /*
         * Gauss-Kronrod-Patterson quadrature coefficients for use in quadpack routine qng.
         *
         * These coefficients were calculated with 101 decimal digit arithmetic by L. W. Fullerton, Bell Labs, Nov 1981.
         */

        // abscissae of the 10-point Gauss-Kronrod-Patterson quadrature formula
        static constexpr const std::array<double, 10> x_10
        {
            -0.973906528517171720077964012084452,
            -0.865063366688984510732096688423493,
            -0.679409568299024406234327365114874,
            -0.433395394129247190799265943165784,
            -0.148874338981631210884826001129720,
            +0.148874338981631210884826001129720,
            +0.433395394129247190799265943165784,
            +0.679409568299024406234327365114874,
            +0.865063366688984510732096688423493,
            +0.973906528517171720077964012084452
        };

        // weights of the 10-point Gauss-Kronrod-Patterson quadrature formula
        static constexpr const std::array<double, 10> w_1_10
        {
            +0.066671344308688137593568809893332,
            +0.149451349150580593145776339657697,
            +0.219086362515982043995534934228163,
            +0.269266719309996355091226921569469,
            +0.295524224714752870173892994651338,
            +0.295524224714752870173892994651338,
            +0.269266719309996355091226921569469,
            +0.219086362515982043995534934228163,
            +0.149451349150580593145776339657697,
            +0.066671344308688137593568809893332
        };

        // weights of the first 10 points of the 21-point Gauss-Kronrod-Patterson quadrature formula
        static constexpr const std::array<double, 10> w_2_10
        {
            +0.032558162307964727478818972459390,
            +0.075039674810919952767043140916190,
            +0.109387158802297641899210590325805,
            +0.134709217311473325928054001771707,
            +0.147739104901338491374841515972068,
            +0.147739104901338491374841515972068,
            +0.134709217311473325928054001771707,
            +0.109387158802297641899210590325805,
            +0.075039674810919952767043140916190,
            +0.032558162307964727478818972459390
        };

        // weights of the first 10 points of the 43-point Gauss-Kronrod-Patterson quadrature formula
        static constexpr const std::array<double, 10> w_3_10
        {
            +0.016296734289666564924281974617663,
            +0.037522876120869501461613795898115,
            +0.054694902058255442147212685465005,
            +0.067355414609478086075553166302174,
            +0.073870199632393953432140695251367,
            +0.073870199632393953432140695251367,
            +0.067355414609478086075553166302174,
            +0.054694902058255442147212685465005,
            +0.037522876120869501461613795898115,
            +0.016296734289666564924281974617663
        };

        // weights of the first 10 points of the 87-point Gauss-Kronrod-Patterson quadrature formula
        static constexpr const std::array<double, 10> w_4_10
        {
            +0.008148377384149172900002878448190,
            +0.018761438201562822243935059003794,
            +0.027347451050052286161582829741283,
            +0.033677707311637930046581056957588,
            +0.036935099820427907614589586742499,
            +0.036935099820427907614589586742499,
            +0.033677707311637930046581056957588,
            +0.027347451050052286161582829741283,
            +0.018761438201562822243935059003794,
            +0.008148377384149172900002878448190
        };

        // abscissae of the 21-point Gauss-Kronrod-Patterson quadrature formula not contained in the 10-point formula
        static constexpr const std::array<double, 11> x_21
        {
            -0.995657163025808080735527280689003,
            -0.930157491355708226001207180059508,
            -0.780817726586416897063717578345042,
            -0.562757134668604683339000099272694,
            -0.294392862701460198131126603103866,
             0.0,
            +0.294392862701460198131126603103866,
            +0.562757134668604683339000099272694,
            +0.780817726586416897063717578345042,
            +0.930157491355708226001207180059508,
            +0.995657163025808080735527280689003
        };

        // weights of the 21-point Gauss-Kronrod-Patterson quadrature formula for abscissae not contained in the 10-point formula
        static constexpr const std::array<double, 11> w_2_21
        {
            +0.011694638867371874278064396062192,
            +0.054755896574351996031381300244580,
            +0.093125454583697605535065465083366,
            +0.123491976262065851077958109831074,
            +0.142775938577060080797094273138717,
            +0.149445554002916905664936468389821,
            +0.142775938577060080797094273138717,
            +0.123491976262065851077958109831074,
            +0.093125454583697605535065465083366,
            +0.054755896574351996031381300244580,
            +0.011694638867371874278064396062192
        };

        // weights of the 43-point Gauss-Kronrod-Patterson quadrature formula for abscissae not contained in the 21-point formula
        static constexpr const std::array<double, 11> w_3_21
        {
            +0.005768556059769796184184327908655,
            +0.027371890593248842081276069289151,
            +0.046560826910428830743339154433824,
            +0.061744995201442564496240336030883,
            +0.071387267268693397768559114425516,
            +0.074722147517403005594425168280423,
            +0.071387267268693397768559114425516,
            +0.061744995201442564496240336030883,
            +0.046560826910428830743339154433824,
            +0.027371890593248842081276069289151,
            +0.005768556059769796184184327908655
        };

        // weights of the 87-point Gauss-Kronrod-Patterson quadrature formula for abscissae not contained in the 43-point formula
        static constexpr const std::array<double, 11> w_4_21
        {
            +0.002884872430211530501334156248695,
            +0.013685946022712701888950035273128,
            +0.023280413502888311123409291030404,
            +0.030872497611713358675466394126442,
            +0.035693633639418770719351355457044,
            +0.037361073762679023410321241766599,
            +0.035693633639418770719351355457044,
            +0.030872497611713358675466394126442,
            +0.023280413502888311123409291030404,
            +0.013685946022712701888950035273128,
            +0.002884872430211530501334156248695,
        };

        // abscissae of the 43-point Gauss-Kronrod-Patterson quadrature formula not contained in the 21-point formula
        static constexpr const std::array<double, 22> x_43
        {
            -0.999333360901932081394099323919911,
            -0.987433402908088869795961478381209,
            -0.954807934814266299257919200290473,
            -0.900148695748328293625099494069092,
            -0.825198314983114150847066732588520,
            -0.732148388989304982612354848755461,
            -0.622847970537725238641159120344323,
            -0.499479574071056499952214885499755,
            -0.364901661346580768043989548502644,
            -0.222254919776601296498260928066212,
            -0.074650617461383322043914435796506,
            +0.074650617461383322043914435796506,
            +0.222254919776601296498260928066212,
            +0.364901661346580768043989548502644,
            +0.499479574071056499952214885499755,
            +0.622847970537725238641159120344323,
            +0.732148388989304982612354848755461,
            +0.825198314983114150847066732588520,
            +0.900148695748328293625099494069092,
            +0.954807934814266299257919200290473,
            +0.987433402908088869795961478381209,
            +0.999333360901932081394099323919911
        };

        // weights of the 43-point Gauss-Kronrod-Patterson quadrature formula for abscissae not contained in the 21-point formula
        static constexpr const std::array<double, 22> w_3_43
        {
            +0.001844477640212414100389106552965,
            +0.010798689585891651740465406741293,
            +0.021895363867795428102523123075149,
            +0.032597463975345689443882222526137,
            +0.042163137935191811847627924327955,
            +0.050741939600184577780189020092084,
            +0.058379395542619248375475369330206,
            +0.064746404951445885544689259517511,
            +0.069566197912356484528633315038405,
            +0.072824441471833208150939535192842,
            +0.074507751014175118273571813842889,
            +0.074507751014175118273571813842889,
            +0.072824441471833208150939535192842,
            +0.069566197912356484528633315038405,
            +0.064746404951445885544689259517511,
            +0.058379395542619248375475369330206,
            +0.050741939600184577780189020092084,
            +0.042163137935191811847627924327955,
            +0.032597463975345689443882222526137,
            +0.021895363867795428102523123075149,
            +0.010798689585891651740465406741293,
            +0.001844477640212414100389106552965
        };

        // weights of the 87-point Gauss-Kronrod-Patterson quadrature formula for abscissae not contained in the 21-point formula
        static constexpr const std::array<double, 22> w_4_43
        {
            +0.000915283345202241360843392549948,
            +0.005399280219300471367738743391053,
            +0.010947679601118931134327826856808,
            +0.016298731696787335262665703223280,
            +0.021081568889203835112433060188190,
            +0.025370969769253827243467999831710,
            +0.029189697756475752501446154084920,
            +0.032373202467202789685788194889595,
            +0.034783098950365142750781997949596,
            +0.036412220731351787562801163687577,
            +0.037253875503047708539592001191226,
            +0.037253875503047708539592001191226,
            +0.036412220731351787562801163687577,
            +0.034783098950365142750781997949596,
            +0.032373202467202789685788194889595,
            +0.029189697756475752501446154084920,
            +0.025370969769253827243467999831710,
            +0.021081568889203835112433060188190,
            +0.016298731696787335262665703223280,
            +0.010947679601118931134327826856808,
            +0.005399280219300471367738743391053,
            +0.000915283345202241360843392549948
        };

        // abscissae of the 87-point Gauss-Kronrod-Patterson quadrature formula not contained in the 43-point formula
        static constexpr const std::array<double, 44> x_87
        {
            -0.999902977262729234490529830591582,
            -0.997989895986678745427496322365960,
            -0.992175497860687222808523352251425,
            -0.981358163572712773571916941623894,
            -0.965057623858384619128284110607926,
            -0.943167613133670596816416634507426,
            -0.915806414685507209591826430720050,
            -0.883221657771316501372117548744163,
            -0.845710748462415666605902011504855,
            -0.803557658035230982788739474980964,
            -0.757005730685495558328942793432020,
            -0.706273209787321819824094274740840,
            -0.651589466501177922534422205016736,
            -0.593223374057961088875273770349144,
            -0.531493605970831932285268948562671,
            -0.466763623042022844871966781659270,
            -0.399424847859218804732101665817923,
            -0.329874877106188288265053371824597,
            -0.258503559202161551802280975429025,
            -0.185695396568346652015917141167606,
            -0.111842213179907468172398359241362,
            -0.037352123394619870814998165437704,
            +0.037352123394619870814998165437704,
            +0.111842213179907468172398359241362,
            +0.185695396568346652015917141167606,
            +0.258503559202161551802280975429025,
            +0.329874877106188288265053371824597,
            +0.399424847859218804732101665817923,
            +0.466763623042022844871966781659270,
            +0.531493605970831932285268948562671,
            +0.593223374057961088875273770349144,
            +0.651589466501177922534422205016736,
            +0.706273209787321819824094274740840,
            +0.757005730685495558328942793432020,
            +0.803557658035230982788739474980964,
            +0.845710748462415666605902011504855,
            +0.883221657771316501372117548744163,
            +0.915806414685507209591826430720050,
            +0.943167613133670596816416634507426,
            +0.965057623858384619128284110607926,
            +0.981358163572712773571916941623894,
            +0.992175497860687222808523352251425,
            +0.997989895986678745427496322365960,
            +0.999902977262729234490529830591582
        };

        // weights of the 87-point Gauss-Kronrod-Patterson quadrature formula for abscissae not contained in the 43-point formula
        static constexpr const std::array<double, 44> w_4_87
        {
            0.000274145563762072350016527092881,
            0.001807124155057942948341311753254,
            0.004096869282759164864458070683480,
            0.006758290051847378699816577897424,
            0.009549957672201646536053581325377,
            0.012329447652244853694626639963780,
            0.015010447346388952376697286041943,
            0.017548967986243191099665352925900,
            0.019938037786440888202278192730714,
            0.022194935961012286796332102959499,
            0.024339147126000805470360647041454,
            0.026374505414839207241503786552615,
            0.028286910788771200659968002987960,
            0.030052581128092695322521110347341,
            0.031646751371439929404586051078883,
            0.033050413419978503290785944862689,
            0.034255099704226061787082821046821,
            0.035262412660156681033782717998428,
            0.036076989622888701185500318003895,
            0.036698604498456094498018047441094,
            0.037120549269832576114119958413599,
            0.037334228751935040321235449094698,
            0.037334228751935040321235449094698,
            0.037120549269832576114119958413599,
            0.036698604498456094498018047441094,
            0.036076989622888701185500318003895,
            0.035262412660156681033782717998428,
            0.034255099704226061787082821046821,
            0.033050413419978503290785944862689,
            0.031646751371439929404586051078883,
            0.030052581128092695322521110347341,
            0.028286910788771200659968002987960,
            0.026374505414839207241503786552615,
            0.024339147126000805470360647041454,
            0.022194935961012286796332102959499,
            0.019938037786440888202278192730714,
            0.017548967986243191099665352925900,
            0.015010447346388952376697286041943,
            0.012329447652244853694626639963780,
            0.009549957672201646536053581325377,
            0.006758290051847378699816577897424,
            0.004096869282759164864458070683480,
            0.001807124155057942948341311753254,
            0.000274145563762072350016527092881,
        };
    }

    template <unsigned dim_>
    std::array<double, dim_>
    integrate(const std::function<std::array<double, dim_> (const double &)> & f,
              const double & a, const double & b,
              const custom::Config & config)
    {
        using value_type = std::array<double, dim_>;

        const double jacobian = (b - a) / 2.0;
        const double x_zero = (b + a) / 2.0;

        std::array<value_type, 10> f_10;
        for (unsigned i = 0 ; i < 10 ; ++i)
        {
            f_10[i] = f(custom::x_10[i] * jacobian + x_zero);
        }

        std::array<value_type, 11> f_21;
        for (unsigned i = 0 ; i < 11 ; ++i)
        {
            f_21[i] = f(custom::x_21[i] * jacobian + x_zero);
        }

        auto accumulator = [](const value_type & a, const value_type & b) -> value_type
        {
            return a + b;
        };

        auto multiplication = [](const double & a, const value_type & b) -> value_type
        {
            return a * b;
        };

        value_type result_1; // result of the 10-point Gauss-Kronrod-Patterson quadrature formula
        result_1.fill(0.0);
        result_1 = std::inner_product(custom::w_1_10.begin(), custom::w_1_10.end(), f_10.begin(), result_1, accumulator, multiplication);
        result_1 = jacobian * result_1;

        value_type result_2; // result of the 21-point Gauss-Kronrod-Patterson quadrature formula
        result_2.fill(0.0);
        result_2 = std::inner_product(custom::w_2_10.begin(), custom::w_2_10.end(), f_10.begin(), result_2, accumulator, multiplication);
        result_2 = std::inner_product(custom::w_2_21.begin(), custom::w_2_21.end(), f_21.begin(), result_2, accumulator, multiplication);
        result_2 = jacobian * result_2;

        auto absolute_accumulator = [](const value_type & a, const value_type & b) -> value_type
        {
            value_type retval;
            for (unsigned i = 0 ; i < dim_ ; ++i)
            {
                retval[i] = std::fabs(a[i]) + std::fabs(b[i]);
            }

            return retval;
        };

        std::array<value_type, 10> fmq2_10; // f(x) - result_2 / (b - a)
        for (unsigned i = 0 ; i < 10 ; ++i)
        {
            fmq2_10[i] = f_10[i] - (1.0 / (b - a)) * result_2;
        }
        std::array<value_type, 11> fmq2_21;
        for (unsigned i = 0 ; i < 11 ; ++i)
        {
            fmq2_21[i] = f_21[i] - (1.0 / (b - a)) * result_2;
        }
        value_type result_absfmq_2; // result of the 21-point Gauss-Kronrod-Patterson quadrature formula for abs(f - result_2 / (b - a))
        result_absfmq_2.fill(0.0);
        result_absfmq_2 = std::inner_product(custom::w_2_10.begin(), custom::w_2_10.end(), fmq2_10.begin(), result_absfmq_2, absolute_accumulator, multiplication);
        result_absfmq_2 = std::inner_product(custom::w_2_21.begin(), custom::w_2_21.end(), fmq2_21.begin(), result_absfmq_2, absolute_accumulator, multiplication);
        result_absfmq_2 = jacobian * result_absfmq_2;

        value_type error_2 = result_2 - result_1;

        bool return_result_2 = true;
        for (unsigned i = 0 ; i < dim_ ; ++i)
        {
            double error = std::fabs(error_2[i]);

            if (error <= std::numeric_limits<double>::min())
            {
                continue;
            }

            if ((result_absfmq_2[i] > 0) && (error > 0))
            {
                error = result_absfmq_2[i] * std::min(1.0, std::pow(200.0 * error / result_absfmq_2[i], 1.5));
            }

            if ((error >= config.epsabs()) && (error >= std::fabs(config.epsrel() * result_2[i])))
            {
                std::cerr << "error_2[" << i << "] = " << error << ">= max(" << config.epsabs() << "," << std::fabs(result_2[i] * config.epsrel()) << ")" << std::endl;
                std::cerr << "result_2[" << i << "] = " << result_2[i] << std::endl;
                return_result_2 = false;
            }
        }

        if (return_result_2)
        {
            return result_2;
        }

        std::array<value_type, 22> f_43;
        for (unsigned i = 0 ; i < 22 ; ++i)
        {
            f_43[i] = f(custom::x_43[i] * jacobian + x_zero);
        }

        value_type result_3; // result of the 43-point Gauss-Kronrod-Patterson quadrature formula
        result_3.fill(0.0);
        result_3 = std::inner_product(custom::w_3_10.begin(), custom::w_3_10.end(), f_10.begin(), result_3, accumulator, multiplication);
        result_3 = std::inner_product(custom::w_3_21.begin(), custom::w_3_21.end(), f_21.begin(), result_3, accumulator, multiplication);
        result_3 = std::inner_product(custom::w_3_43.begin(), custom::w_3_43.end(), f_43.begin(), result_3, accumulator, multiplication);
        result_3 = jacobian * result_3;

        std::array<value_type, 10> fmq3_10;
        for (unsigned i = 0 ; i < 10 ; ++i)
        {
            fmq3_10[i] = f_10[i] - (1.0 / (b - a)) * result_3;
        }
        std::array<value_type, 11> fmq3_21;
        for (unsigned i = 0 ; i < 11 ; ++i)
        {
            fmq3_21[i] = f_21[i] - (1.0 / (b - a)) * result_3;
        }
        std::array<value_type, 22> fmq3_43;
        for (unsigned i = 0 ; i < 22 ; ++i)
        {
            fmq3_43[i] = f_43[i] - (1.0 / (b - a)) * result_3;
        }

        value_type result_absfmq_3; // result of the 43-point Gauss-Kronrod-Patterson quadrature formula for abs(f - result_3 / (b - a))
        result_absfmq_3.fill(0.0);
        result_absfmq_3 = std::inner_product(custom::w_3_10.begin(), custom::w_3_10.end(), fmq3_10.begin(), result_absfmq_3, absolute_accumulator, multiplication);
        result_absfmq_3 = std::inner_product(custom::w_3_21.begin(), custom::w_3_21.end(), fmq3_21.begin(), result_absfmq_3, absolute_accumulator, multiplication);
        result_absfmq_3 = std::inner_product(custom::w_3_43.begin(), custom::w_3_43.end(), fmq3_43.begin(), result_absfmq_3, absolute_accumulator, multiplication);
        result_absfmq_3 = jacobian * result_absfmq_3;

        value_type error_3 = result_3 - result_2;

        bool return_result_3 = true;
        for (unsigned i = 0 ; i < dim_ ; ++i)
        {
            double error = std::fabs(error_3[i]);

            if (error <= std::numeric_limits<double>::min())
            {
                continue;
            }

            if ((result_absfmq_3[i] > 0) && (error > 0))
            {
                error = result_absfmq_3[i] * std::min(1.0, std::pow(200.0 * error / result_absfmq_3[i], 1.5));
            }

            if ((error >= config.epsabs()) && (error >= std::fabs(config.epsrel() * result_3[i])))
            {
                std::cerr << "error_3[" << i << "] = " << error << ">= max(" << config.epsabs() << "," << std::fabs(result_3[i] * config.epsrel()) << ")" << std::endl;
                std::cerr << "result_3[" << i << "] = " << result_3[i] << std::endl;
                return_result_3 = false;
            }
        }

        if (return_result_3)
        {
            return result_3;
        }

        std::array<value_type, 44> f_87;
        for (unsigned i = 0 ; i < 44 ; ++i)
        {
            f_87[i] = f(custom::x_87[i] * jacobian + x_zero);
        }

        value_type result_4; // result of the 87-point Gauss-Kronrod-Patterson quadrature formula
        result_4.fill(0.0);
        result_4 = std::inner_product(custom::w_4_10.begin(), custom::w_4_10.end(), f_10.begin(), result_4, accumulator, multiplication);
        result_4 = std::inner_product(custom::w_4_21.begin(), custom::w_4_21.end(), f_21.begin(), result_4, accumulator, multiplication);
        result_4 = std::inner_product(custom::w_4_43.begin(), custom::w_4_43.end(), f_43.begin(), result_4, accumulator, multiplication);
        result_4 = std::inner_product(custom::w_4_87.begin(), custom::w_4_87.end(), f_87.begin(), result_4, accumulator, multiplication);
        result_4 = jacobian * result_4;

        std::array<value_type, 10> fmq4_10;
        for (unsigned i = 0 ; i < 10 ; ++i)
        {
            fmq4_10[i] = f_10[i] - (1.0 / (b - a)) * result_4;
        }
        std::array<value_type, 11> fmq4_21;
        for (unsigned i = 0 ; i < 11 ; ++i)
        {
            fmq4_21[i] = f_21[i] - (1.0 / (b - a)) * result_4;
        }
        std::array<value_type, 22> fmq4_43;
        for (unsigned i = 0 ; i < 22 ; ++i)
        {
            fmq4_43[i] = f_43[i] - (1.0 / (b - a)) * result_4;
        }
        std::array<value_type, 44> fmq4_87;
        for (unsigned i = 0 ; i < 44 ; ++i)
        {
            fmq4_87[i] = f_87[i] - (1.0 / (b - a)) * result_4;
        }

        value_type result_absfmq_4; // result of the 87-point Gauss-Kronrod-Patterson quadrature formula for abs(f - result_4 / (b - a))
        result_absfmq_4.fill(0.0);
        result_absfmq_4 = std::inner_product(custom::w_4_10.begin(), custom::w_4_10.end(), fmq4_10.begin(), result_absfmq_4, absolute_accumulator, multiplication);
        result_absfmq_4 = std::inner_product(custom::w_4_21.begin(), custom::w_4_21.end(), fmq4_21.begin(), result_absfmq_4, absolute_accumulator, multiplication);
        result_absfmq_4 = std::inner_product(custom::w_4_43.begin(), custom::w_4_43.end(), fmq4_43.begin(), result_absfmq_4, absolute_accumulator, multiplication);
        result_absfmq_4 = std::inner_product(custom::w_4_87.begin(), custom::w_4_87.end(), fmq4_87.begin(), result_absfmq_4, absolute_accumulator, multiplication);
        result_absfmq_4 = jacobian * result_absfmq_4;

        value_type error_4 = result_4 - result_3;

        bool return_result_4 = true;
        for (unsigned i = 0 ; i < dim_ ; ++i)
        {
            double error = std::fabs(error_4[i]);

            if (error <= std::numeric_limits<double>::min())
            {
                continue;
            }

            if ((result_absfmq_4[i] > 0) && (error > 0))
            {
                error = result_absfmq_4[i] * std::min(1.0, std::pow(200.0 * error / result_absfmq_4[i], 1.5));
            }

            if ((error >= config.epsabs()) && (error >= std::fabs(config.epsrel() * result_4[i])))
            {
                std::cerr << "error_4[" << i << "] = " << error << ">= max(" << config.epsabs() << "," << std::fabs(result_4[i] * config.epsrel()) << ")" << std::endl;
                std::cerr << "result_4[" << i << "] = " << result_4[i] << std::endl;
                return_result_4 = false;
            }
        }

        if (return_result_4)
        {
            return result_4;
        }

        //throw IntegrationError("Gauss-Kronrod-Patterson 87-point quadrature failed to converge");
        return result_4; // to silence compiler warning
    }

    template <std::size_t k> std::array<double, k> integrate1D(const std::function<std::array<double, k> (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        if (n < 16)
            n = 16;

        // step width
        double h = (b - a) / n;

        // evaluate function for every sampling point
        std::vector<std::array<double, k>> y;
        for (unsigned i = 0 ; i < n + 1 ; ++i)
        {
            y.push_back(f(a + i * h));
        }

        std::array<double, k> Q0; Q0.fill(0.0);
        std::array<double, k> Q1; Q1.fill(0.0);
        std::array<double, k> Q2; Q2.fill(0.0);

        for (unsigned i = 0 ; i < n / 8 ; ++i)
        {
            Q0 = Q0 + y[8 * i] + 4.0 * y[8 * i + 4] + y[8 * i + 4];
        }
        for (unsigned i = 0 ; i < n / 4 ; ++i)
        {
            Q1 = Q1 + y[4 * i] + 4.0 * y[4 * i + 2] + y[4 * i + 4];
        }
        for (unsigned i = 0 ; i < n / 2 ; ++i)
        {
            Q2 = Q2 + y[2 * i] + 4.0 * y[2 * i + 1] + y[2 * i + 2];
        }

        Q0 = (h / 3.0 * 4.0) * Q0;
        Q1 = (h / 3.0 * 2.0) * Q1;
        Q2 = (h / 3.0) * Q2;

        std::array<double, k> denom = Q0 + Q2 - 2.0 * Q1;
        std::array<double, k> num = Q2 - Q1;
        std::array<double, k> correction = divide(mult(num, num), denom);

        bool correction_valid = true;
        for (unsigned i = 0 ; i < k ; ++i)
        {
            if (std::isnan(correction[i]))
            {
                correction_valid = false;
                break;
            }
        }

        if (!correction_valid)
        {
            return Q2;
        }
        else
        {
            bool correction_small = true;

            for (unsigned i = 0 ; i < k ; ++i)
            {
                if ((abs(correction[i] / Q2[i])) > 1.0)
                {
                    correction_small = false;
                    break;
                }
            }

            if (correction_small)
            {
                return Q2 - correction;
            }
            else
            {
#if 0
                std::cerr << "Q0 = " << Q0 << std::endl;
                std::cerr << "Q1 = " << Q1 << std::endl;
                std::cerr << "Q2 = " << Q2 << std::endl;
                std::cerr << "Reintegrating with twice the number of data points" << std::endl;
#endif
                return integrate1D(f, 2 * n, a, b);
            }
        }
    }

    namespace cubature
    {

        template <size_t dim_>
        int scalar_integrand(unsigned ndim , const double *x, void *data,
                      unsigned fdim , double *fval)
        {
            assert(ndim == dim_);
            assert(fdim == 1);

            auto& f = *static_cast<cubature::fdd<dim_> *>(data);
            // TODO use std::array_view once available
            std::array<double, dim_> args;
            std::copy(x, x + dim_, args.data());
            *fval = f(args);

            return 0;
        }

    }

    template <size_t dim_>
    double integrate(const cubature::fdd<dim_> & f,
                     const std::array<double, dim_> &a,
                     const std::array<double, dim_> &b,
                     const cubature::Config &config)
    {
        // TODO Support infinite intervals by param trafo? Not for now.
        constexpr unsigned nintegrands = 1;
        double res;
        double err;
        if (hcubature(nintegrands, &cubature::scalar_integrand<dim_>,
                      &const_cast<cubature::fdd<dim_>&>(f), dim_, a.data(), b.data(),
                      config.maxeval(), config.epsabs(), config.epsrel(), ERROR_L2, &res, &err))
        {
            throw IntegrationError("hcubature failed");
        }

        return res;
    }

}

#endif
