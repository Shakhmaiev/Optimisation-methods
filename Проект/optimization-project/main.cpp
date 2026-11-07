
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "src/variant.h"
#include "src/image.h"
#include "src/metrics.h"
#include "src/filters.h"
#include "src/oracle.h"
#include "src/lab1_svenn.h"
#include "src/lab2_search.h"
#include "src/lab3_objective.h"
#include "src/lab3_nelder_mead.h"
#include "src/lab3_powell.h"
#include "src/lab4_tests.h"
#include "src/lab4_lambda.h"
#include "src/lab4_reference.h"
#include "src/lab4_gradient.h"
#include "src/lab5_constrained.h"
using namespace std;

int main() {

    cout << fixed << setprecision(8);

    cout << "================================\n";
    cout << "     OPTIMIZATION PROJECT\n";
    cout << "          VARIANT 15\n";
    cout << "================================\n";

    cout << "1 - Lab 2: Main results\n";
    cout << "2 - Lab 3: Nelder Mead\n";
    cout << "3 - Lab 3: Powell\n";
    cout << "4 - Lab 3: Experiment 1\n";
    cout << "5 - Lab 3: Experiment 3\n";
    cout << "6 - Lab 3: Experiment 4\n";
    cout << "7 - Lab 3: Contour data\n";
    cout << "8 - Lab 4: Operators test\n";
    cout << "9 - Lab 4: Lambda test\n";
    cout << "10 - Lab 4: Lambda optimization\n";
    cout << "11 - Lab 4: Reference solution\n";
    cout << "12 - Lab 4: Gradient methods\n";
    cout << "13 - Lab 4: Experiment 1\n";
    cout << "14 - Lab 4: Experiment 2\n";
    cout << "15 - Lab 4: Experiment 3\n"; 
    cout << "16 - Lab 5: Raw image test\n";
    cout << "17 - Lab 5: Projected gradient\n";
    cout << "18 - Lab 5: Penalty method\n";
    cout << "19 - Lab 5: Experiments\n";
    cout << "20 - Lab 5: FISTA\n";
    cout << "21 - Lab 1: Svenn and step experiments\n";
    cout << "22 - Lab 2: Iterations and epsilon experiments\n";
    cout << "23 - Lab 5: Unconstrained baseline\n";


    
    cout << "0 - Exit\n";

    cout << "\nChoose: ";

    int choice;
    cin >> choice;

    if (choice == 0)
        return 0;

    // Загальні дані
    Image clean = clean_image(SIZE);
    NoisyData data = add_noise(clean, NOISE_SIGMA, SEED);

    auto f = [&](double sigma) {
        double s = clamp(sigma, 0.05, 8.0);
        Image restored = gaussian_blur(data.noisy, s);

        return mse(restored, clean)
             + 1000.0 * (sigma - s) * (sigma - s);
    };

    Oracle oracle(f);

    // Лабораторна 2
    if (choice == 1) {

        SvennResult interval = svenn(
            oracle, X0, DELTA0
        );

        oracle.reset();

        SearchResult d = dichotomy(
            oracle, interval.a, interval.b, EPS
        );

        oracle.reset();

        SearchResult g = golden_section(
            oracle, interval.a, interval.b, EPS
        );

        cout << "\nLAB 2 RESULTS\n";

        cout << "Dichotomy sigma = " << d.x << endl;
        cout << "Dichotomy MSE = " << d.f << endl;
        cout << "Dichotomy calls = " << d.n_f << endl;

        cout << "Golden sigma = " << g.x << endl;
        cout << "Golden MSE = " << g.f << endl;
        cout << "Golden calls = " << g.n_f << endl;

        return 0;
    }

    // Початкові параметри лабораторної 3
    // Використовуємо результат лабораторної 2
    // Derive the spatial-filter starting point from Lab 2,
    // rather than copying a previously observed result.
    Oracle init_oracle(f);
    SvennResult init_interval = svenn(init_oracle, X0, DELTA0);
    init_oracle.reset();
    const double sigma_gaussian = golden_section(
        init_oracle, init_interval.a, init_interval.b, EPS
    ).x;

    double sigma_s = clamp(
        2.0 * sigma_gaussian, 0.5, 4.0
    );

    double sigma_r = clamp(
        0.15 + 0.1 * NOISE_SIGMA,
        0.02, 0.9
    );

    Oracle2 oracle2([&](double s, double r) {
        return bilateral_objective(
            data.noisy, clean, s, r
        );
    });


    // Lambda is recalculated using the practical's independent search settings.
    auto calculate_lambda = [&]() -> double {
        Oracle lam_oracle([&](double lam) {
            return lambda_objective(data.noisy, clean, lam);
        });
        SvennResult bracket = svenn(lam_oracle, LAMBDA0, 0.2, 1e-4, 50.0);
        lam_oracle.reset();
        return golden_section(lam_oracle, bracket.a, bracket.b, 1e-3).x;
    };

    // LAB 3 - NELDER MEAD
    if (choice == 2) {

        NMResult nm = nelder_mead(
            oracle2, sigma_s, sigma_r, SIMPLEX_STEP
        );

        Image restored = bilateral_filter(
            data.noisy, nm.a, nm.b
        );

        cout << "\nLAB 3 - NELDER MEAD\n";
        cout << "Sigma spatial = " << nm.a << endl;
        cout << "Sigma range = " << nm.b << endl;
        cout << "MSE = " << nm.f << endl;
        cout << "Iterations = " << nm.iterations << endl;
        cout << "Calls = " << nm.calls << endl;
        cout << "PSNR = " << psnr(restored, clean) << endl;

        save_pgm(restored, "outputs/lab3_nelder_mead.pgm");


    ofstream nm_csv("outputs/lab3_nelder_mead.csv");

    nm_csv << "k;operation;sigma_s;sigma_r;MSE;size;calls\n";

    for (const auto& row : nm.rows) {
        nm_csv << row.k << ";"
               << row.operation << ";"
               << row.a << ";"
               << row.b << ";"
               << row.f << ";"
               << row.size << ";"
               << row.calls << "\n";
    }

    nm_csv.close();


        return 0;
    }

    // LAB 3 - POWELL
    if (choice == 3) {

        PowellResult pw = powell(
            oracle2, sigma_s, sigma_r
        );

        Image restored = bilateral_filter(
            data.noisy, pw.a, pw.b
        );

        cout << "\nLAB 3 - POWELL\n";
        cout << "Sigma spatial = " << pw.a << endl;
        cout << "Sigma range = " << pw.b << endl;
        cout << "MSE = " << pw.f << endl;
        cout << "Cycles = " << pw.cycles << endl;
        cout << "Calls = " << pw.calls << endl;
        cout << "PSNR = " << psnr(restored, clean) << endl;

        save_pgm(restored, "outputs/lab3_powell.pgm");


    ofstream pw_csv("outputs/lab3_powell.csv");

    pw_csv << "cycle;sigma_s;sigma_r;MSE;shift;calls\n";

    for (const auto& row : pw.rows) {
        pw_csv << row.cycle << ";"
               << row.a << ";"
               << row.b << ";"
               << row.f << ";"
               << row.shift << ";"
               << row.calls << "\n";
    }

    pw_csv.close();


        return 0;
    }

    // EXPERIMENT 1 - SIMPLEX SIZE
    if (choice == 4) {

        cout << "\nLAB 3 - EXPERIMENT 1\n";

        double steps[] = {
            SIMPLEX_STEP / 3.0,
            SIMPLEX_STEP,
            SIMPLEX_STEP * 3.0
        };

        ofstream csv("outputs/lab3_simplex_experiment.csv");
        csv << "step;sigma_s;sigma_r;MSE;iterations;calls\n";

        for (double step : steps) {

            oracle2.reset();

            NMResult result = nelder_mead(
                oracle2, sigma_s, sigma_r, step
            );

            cout << "Step = " << step
                 << " Sigma s = " << result.a
                 << " Sigma r = " << result.b
                 << " MSE = " << result.f
                 << " Iterations = " << result.iterations
                 << " Calls = " << result.calls << endl;

            csv << step << ";"
                << result.a << ";"
                << result.b << ";"
                << result.f << ";"
                << result.iterations << ";"
                << result.calls << "\n";
        }

        return 0;
    }

    // EXPERIMENT 3 - POWELL TOLERANCE
    if (choice == 5) {

        cout << "\nLAB 3 - EXPERIMENT 3\n";

        double tolerances[] = {1e-4, 1e-7};

        ofstream csv("outputs/lab3_tolerance_experiment.csv");
        csv << "tol;sigma_s;sigma_r;MSE;cycles;calls\n";

        for (double tol : tolerances) {

            oracle2.reset();

            PowellResult result = powell(
                oracle2, sigma_s, sigma_r, tol
            );

            cout << "Tolerance = " << tol
                 << " Sigma s = " << result.a
                 << " Sigma r = " << result.b
                 << " MSE = " << result.f
                 << " Cycles = " << result.cycles
                 << " Calls = " << result.calls << endl;

            csv << tol << ";"
                << result.a << ";"
                << result.b << ";"
                << result.f << ";"
                << result.cycles << ";"
                << result.calls << "\n";
        }

        return 0;
    }

    // EXPERIMENT 4 - ONE VS TWO PARAMETERS
    if (choice == 6) {

        cout << "\nLAB 3 - EXPERIMENT 4\n";

        double fixed_s = sigma_gaussian;

        Oracle range_oracle([&](double r) {
            return bilateral_objective(
                data.noisy, clean, fixed_s, r
            );
        });

        SearchResult one = golden_section(
            range_oracle, 0.01, 1.0, EPS
        );

        Image one_image = bilateral_filter(
            data.noisy, fixed_s, one.x
        );

        // Результат двовимірної оптимізації,
        // який ми вже отримали
        // Точка зліва від розриву, близька до результату Нелдера–Міда
Oracle2 best_oracle([&](double s, double r) {
    return bilateral_objective(data.noisy, clean, s, r);
});
NMResult best_nm = nelder_mead(best_oracle, sigma_s, sigma_r, SIMPLEX_STEP);
const double best_s = best_nm.a;
const double best_r = best_nm.b;
cout << "\nBOUNDARY CHECK:\n";

for (double s : {1.749999, 1.7499999, 1.75, 1.7500001}) {
    cout << "Sigma s = " << setprecision(10) << s
         << " MSE = "
         << bilateral_objective(data.noisy, clean, s, best_r)
         << endl;
}

        Image two_image = bilateral_filter(
            data.noisy, best_s, best_r
        );

        double mse_two = mse(two_image, clean);
        double psnr_one = psnr(one_image, clean);
        double psnr_two = psnr(two_image, clean);

        cout << "ONE PARAMETER:\n";
        cout << "Sigma s = " << fixed_s << endl;
        cout << "Sigma r = " << one.x << endl;
        cout << "MSE = " << one.f << endl;
        cout << "PSNR = " << psnr_one << endl;

        cout << "\nTWO PARAMETERS:\n";
        cout << "Sigma s = " << best_s << endl;
        cout << "Sigma r = " << best_r << endl;
        cout << "MSE = " << mse_two << endl;
        cout << "PSNR = " << psnr_two << endl;

        cout << "\nMSE difference = "
             << one.f - mse_two << endl;

        cout << "PSNR difference = "
             << psnr_two - psnr_one << endl;

        save_pgm(one_image, "outputs/lab3_one_parameter.pgm");

        ofstream csv("outputs/lab3_one_vs_two.csv");
        csv << "method;sigma_s;sigma_r;MSE;PSNR\n";

        csv << "One parameter;"
            << fixed_s << ";" << one.x << ";"
            << one.f << ";" << psnr_one << "\n";

        csv << "Two parameters;"
            << best_s << ";" << best_r << ";"
            << mse_two << ";" << psnr_two << "\n";

        return 0;
    }


    // LAB 3 - CONTOUR DATA
    // Сітка значень цільової функції

    if (choice == 7) {

        cout << "\nLAB 3 - CONTOUR DATA\n";

        ofstream csv("outputs/lab3_contours.csv");

        csv << fixed << setprecision(12);
        csv << "sigma_s;sigma_r;MSE\n";

        const int N = 40;

        const double s_min = 0.3;
        const double s_max = 5.0;

        const double r_min = 0.01;
        const double r_max = 1.0;

        double best_mse = 1e100;
        double best_s = 0.0;
        double best_r = 0.0;

        for (int i = 0; i < N; i++) {

            double s = s_min +
                (s_max - s_min) * i / (N - 1);

            for (int j = 0; j < N; j++) {

                double r = r_min +
                    (r_max - r_min) * j / (N - 1);

                // peek не збільшує лічильник
                double error = oracle2.peek(s, r);

                csv << s << ";"
                    << r << ";"
                    << error << "\n";

                if (error < best_mse) {
                    best_mse = error;
                    best_s = s;
                    best_r = r;
                }
            }

            cout << "Grid row " << i + 1
                 << "/" << N << endl;
        }

        csv.close();

        cout << "\nGrid completed!\n";
        cout << "Grid points = " << N * N << endl;
        cout << "Best grid sigma s = " << best_s << endl;
        cout << "Best grid sigma r = " << best_r << endl;
        cout << "Best grid MSE = " << best_mse << endl;

        cout << "Saved: outputs/lab3_contours.csv\n";

        return 0;
    }

if (choice == 8) {
    lab4_tests();
    return 0;
}



    // LAB 4 - LAMBDA TEST
    if (choice == 9) {

        cout << "\nLAB 4 - LAMBDA TEST\n";

        cout << "Initial lambda = " << LAMBDA0 << endl;

        double L = 1.0 + 8.0 * LAMBDA0;
        double mu = 1.0;
        double kappa = L / mu;
        double h = 1.0 / L;

        cout << "L = " << L << endl;
        cout << "Mu = " << mu << endl;
        cout << "Kappa = " << kappa << endl;
        cout << "Step h = " << h << endl;
        cout << "Stability limit = " << 2.0 / L << endl;

        double lambdas[] = {
            0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0
        };

        ofstream csv("outputs/lab4_lambda_test.csv");
        csv << fixed << setprecision(12);
        csv << "lambda;MSE;PSNR\n";

        double best_lambda = 0.0;
        double best_mse = 1e100;

        for (double lambda : lambdas) {

            Image restored = tikhonov_solve(
                data.noisy, lambda
            );

            double error = mse(restored, clean);
            double quality = psnr(restored, clean);

            cout << "Lambda = " << lambda
                 << " MSE = " << error
                 << " PSNR = " << quality
                 << endl;

            csv << lambda << ";"
                << error << ";"
                << quality << "\n";

            if (error < best_mse) {
                best_mse = error;
                best_lambda = lambda;
            }
        }

        csv.close();

        cout << "\nBest tested lambda = "
             << best_lambda << endl;

        cout << "Best tested MSE = "
             << best_mse << endl;

        cout << "Saved: outputs/lab4_lambda_test.csv\n";

        return 0;
    }


    // LAB 4 - LAMBDA OPTIMIZATION
    if (choice == 10) {

        cout << "\nLAB 4 - LAMBDA OPTIMIZATION\n";

        // Зовнішня цільова функція
        Oracle lambda_oracle([&](double lambda) {
            return lambda_objective(
                data.noisy, clean, lambda
            );
        });

        // Пошук інтервалу методом Свенна
        SvennResult interval = svenn(
            lambda_oracle,
            LAMBDA0,
            0.2, 0.0001, 50.0
        );

        cout << "\nSVENN:\n";
        cout << "a = " << interval.a << endl;
        cout << "b = " << interval.b << endl;
        cout << "Calls = " << interval.n_f << endl;

        // Уточнення золотим перерізом
        lambda_oracle.reset();

        SearchResult result = golden_section(
            lambda_oracle,
            interval.a,
            interval.b,
            1e-3
        );

        double lambda_star = result.x;

        // Відновлене зображення
        Image restored = tikhonov_solve(
            data.noisy,
            lambda_star
        );

        // Константи оптимізації
        double L = 1.0 + 8.0 * lambda_star;
        double mu = 1.0;
        double kappa = L / mu;
        double h = 1.0 / L;

        cout << "\nGOLDEN SECTION:\n";
        cout << "Optimal lambda = " << lambda_star << endl;
        cout << "MSE = " << mse(restored, clean) << endl;
        cout << "PSNR = " << psnr(restored, clean) << endl;
        cout << "Iterations = " << result.iterations << endl;
        cout << "Calls = " << result.n_f << endl;

        cout << "\nCONSTANTS:\n";
        cout << "L = " << L << endl;
        cout << "Mu = " << mu << endl;
        cout << "Kappa = " << kappa << endl;
        cout << "Step h = " << h << endl;
        cout << "Stability limit = " << 2.0 / L << endl;

        cout << "E* = " << energy(restored, data.noisy, lambda_star) << endl;

        save_pgm(
            restored,
            "outputs/lab4_lambda_optimal.pgm"
        );

        // Таблиця ітерацій
        ofstream csv("outputs/lab4_lambda_search.csv");

        csv << fixed << setprecision(12);
        csv << "k;a;b;x1;x2;f1;f2;length\n";

        for (const auto& row : result.rows) {
            csv << row.k << ";"
                << row.a << ";"
                << row.b << ";"
                << row.x1 << ";"
                << row.x2 << ";"
                << row.f1 << ";"
                << row.f2 << ";"
                << row.length << "\n";
        }

        csv.close();

        cout << "\nSaved: outputs/lab4_lambda_search.csv\n";
        cout << "Saved: outputs/lab4_lambda_optimal.pgm\n";

        return 0;
    }


    // LAB 4 - REFERENCE SOLUTION
    if (choice == 11) {

        cout << "\nLAB 4 - REFERENCE SOLUTION\n";

        const double lambda_star = calculate_lambda();

        ReferenceResult ref = find_reference(
            data.noisy,
            lambda_star
        );

        cout << "Lambda = " << lambda_star << endl;
        cout << "Iterations = " << ref.iterations << endl;
        cout << "Energy E* = " << ref.E << endl;
        cout << "Gradient norm = " << ref.grad_norm << endl;

        cout << "Converged = "
             << (ref.converged ? "YES" : "NO")
             << endl;

        cout << "MSE = " << mse(ref.image, clean) << endl;
        cout << "PSNR = " << psnr(ref.image, clean) << endl;

        save_pgm(
            ref.image,
            "outputs/lab4_reference.pgm"
        );

        cout << "Saved: outputs/lab4_reference.pgm\n";

        return 0;
    }


    // LAB 4 - GRADIENT METHODS
    if (choice == 12) {

        cout << "\nLAB 4 - GRADIENT METHODS\n";

        const double lambda = calculate_lambda();

        double L = 1.0 + 8.0 * lambda;
        double kappa = L;

        // Коефіцієнти інерції
        double beta_hb = pow(
            (sqrt(kappa) - 1.0) /
            (sqrt(kappa) + 1.0), 2.0
        );

        double beta_nag =
            (sqrt(kappa) - 1.0) /
            (sqrt(kappa) + 1.0);

        struct MethodConfig {
            string name;
            string kind;
            double beta;
        };

        vector<MethodConfig> methods = {
            {"Gradient Descent", "gd", 0.0},
            {"Armijo", "armijo", 0.0},
            {"Heavy Ball Optimal", "hb", beta_hb},
            {"Heavy Ball 0.9", "hb", 0.9},
            {"Nesterov Convex", "nag", -1.0},
            {"Nesterov Strong", "nag", beta_nag}
        };

        ofstream summary("outputs/lab4_methods.csv");

        summary << fixed << setprecision(12);
        summary << "method;iterations;calls;grad_calls;"
                   "E;grad_norm;converged;PSNR\n";

        for (const auto& method : methods) {

            GradResult result = gradient_method(
                data.noisy,
                lambda,
                method.kind,
                method.beta
            );

            double quality = psnr(result.image, clean);

            cout << "\nMethod: " << method.name << endl;
            cout << "Iterations = " << result.iterations << endl;
            cout << "Energy = " << result.E << endl;
            cout << "Gradient norm = " << result.grad_norm << endl;
            cout << "Function calls = " << result.calls << endl;
            cout << "Gradient calls = " << result.grad_calls << endl;
            cout << "PSNR = " << quality << endl;
            cout << "Converged = "
                 << (result.converged ? "YES" : "NO")
                 << endl;

            summary << method.name << ";"
                    << result.iterations << ";"
                    << result.calls << ";"
                    << result.grad_calls << ";"
                    << result.E << ";"
                    << result.grad_norm << ";"
                    << (result.converged ? 1 : 0) << ";"
                    << quality << "\n";
        }

        summary.close();

        cout << "\nSaved: outputs/lab4_methods.csv\n";

        return 0;
    }




    // LAB 4 - EXPERIMENT 1
    if (choice == 13) {

        cout << "\nLAB 4 - EXPERIMENT 1\n";

        const double lambda = 20.0;

        const double L = 1.0 + 8.0 * lambda;
        const double mu = 1.0;
        const double kappa = L / mu;

        const double beta_hb = pow(
            (sqrt(L) - sqrt(mu)) /
            (sqrt(L) + sqrt(mu)), 2.0
        );

        const double beta_nag =
            (sqrt(L) - sqrt(mu)) /
            (sqrt(L) + sqrt(mu));

        cout << "Lambda = " << lambda << endl;
        cout << "L = " << L << endl;
        cout << "Kappa = " << kappa << endl;

        struct Config {
            string name;
            string kind;
            double beta;
        };

        vector<Config> methods = {
            {"Gradient Descent", "gd", 0.0},
            {"Heavy Ball Optimal", "hb", beta_hb},
            {"Nesterov Strong", "nag", beta_nag}
        };

        ofstream csv("outputs/lab4_experiment1.csv");
        csv << fixed << setprecision(12);
        csv << "method;iterations;calls;grad_calls;"
               "E;grad_norm;converged\n";

        for (const auto& method : methods) {

            GradResult result = gradient_method(
                data.noisy,
                lambda,
                method.kind,
                method.beta,
                3000,
                1e-6
            );

            cout << "\nMethod: " << method.name << endl;
            cout << "Iterations = " << result.iterations << endl;
            cout << "Energy = " << result.E << endl;
            cout << "Gradient norm = " << result.grad_norm << endl;
            cout << "Function calls = " << result.calls << endl;
            cout << "Converged = "
                 << (result.converged ? "YES" : "NO")
                 << endl;

            csv << method.name << ";"
                << result.iterations << ";"
                << result.calls << ";"
                << result.grad_calls << ";"
                << result.E << ";"
                << result.grad_norm << ";"
                << (result.converged ? 1 : 0)
                << "\n";
        }

        csv.close();

        cout << "\nSaved: outputs/lab4_experiment1.csv\n";

        return 0;
    }



    // LAB 4 - EXPERIMENT 2
    if (choice == 14) {

        cout << "\nLAB 4 - EXPERIMENT 2\n";

        const double lambda = calculate_lambda();
        const double L = 1.0 + 8.0 * lambda;

        const double h_stable = 1.0 / L;
        const double h_unstable = 2.2 / L;

        cout << "Lambda = " << lambda << endl;
        cout << "L = " << L << endl;
        cout << "Stability limit = " << 2.0 / L << endl;
        cout << "Stable step = " << h_stable << endl;
        cout << "Unstable step = " << h_unstable << endl;

        GradResult stable = gradient_method(
            data.noisy,
            lambda,
            "gd",
            0.0,
            50,
            0.0,
            h_stable
        );

        GradResult unstable = gradient_method(
            data.noisy,
            lambda,
            "gd",
            0.0,
            50,
            0.0,
            h_unstable
        );

        ofstream csv("outputs/lab4_experiment2.csv");
        csv << fixed << setprecision(12);
        csv << "k;E_stable;E_unstable\n";

        cout << "\nITERATION COMPARISON:\n";

        const size_t n = min(
            stable.rows.size(),
            unstable.rows.size()
        );

        for (size_t i = 0; i < n; i++) {

            const auto& a = stable.rows[i];
            const auto& b = unstable.rows[i];

            csv << a.k << ";"
                << a.E << ";"
                << b.E << "\n";

            if (a.k % 10 == 0) {
                cout << "k = " << a.k
                     << " E stable = " << a.E
                     << " E unstable = " << b.E
                     << endl;
            }
        }

        // Додаємо значення після 50-го кроку
        csv << stable.iterations << ";"
            << stable.E << ";"
            << unstable.E << "\n";

        cout << "\nFINAL RESULTS:\n";
        cout << "Stable E = " << stable.E << endl;
        cout << "Unstable E = " << unstable.E << endl;

        csv.close();

        cout << "\nSaved: outputs/lab4_experiment2.csv\n";

        return 0;
    }



    // LAB 4 - EXPERIMENT 3
    if (choice == 15) {

        cout << "\nLAB 4 - EXPERIMENT 3\n";

        const double lambda_star = calculate_lambda();
        const int N = 101;

        const double lambda_min = 0.05;
        const double lambda_max = 1.25;

        ofstream csv("outputs/lab4_experiment3.csv");

        csv << fixed << setprecision(12);
        csv << "lambda;MSE;PSNR\n";

        double best_lambda = 0.0;
        double best_mse = 1e100;
        double best_psnr = 0.0;

        for (int i = 0; i < N; i++) {

            double lambda = lambda_min
                + (lambda_max - lambda_min)
                * i / (N - 1);

            Image restored = tikhonov_solve(
                data.noisy, lambda, 300
            );

            double error = mse(restored, clean);
            double quality = psnr(restored, clean);

            csv << lambda << ";"
                << error << ";"
                << quality << "\n";

            if (error < best_mse) {
                best_mse = error;
                best_lambda = lambda;
                best_psnr = quality;
            }
        }

        // Точний результат золотого перерізу
        Image optimal = tikhonov_solve(
            data.noisy, lambda_star, 300
        );

        double optimal_mse = mse(optimal, clean);
        double optimal_psnr = psnr(optimal, clean);

        csv << lambda_star << ";"
            << optimal_mse << ";"
            << optimal_psnr << "\n";

        csv.close();

        cout << "Grid points = " << N << endl;

        cout << "\nBEST GRID RESULT:\n";
        cout << "Lambda = " << best_lambda << endl;
        cout << "MSE = " << best_mse << endl;
        cout << "PSNR = " << best_psnr << endl;

        cout << "\nGOLDEN SECTION RESULT:\n";
        cout << "Lambda = " << lambda_star << endl;
        cout << "MSE = " << optimal_mse << endl;
        cout << "PSNR = " << optimal_psnr << endl;

        cout << "\nSaved: outputs/lab4_experiment3.csv\n";

        return 0;
    }



    // LAB 5 - RAW IMAGE TEST
    if (choice == 16) {

        cout << "\nLAB 5 - RAW IMAGE TEST\n";

        const Image& raw = data.raw;

        Violation v = violation(raw);

        auto limits = minmax_element(
            raw.px.begin(),
            raw.px.end()
        );

        cout << "Pixels = " << raw.px.size() << endl;
        cout << "Outside [0,1] = " << v.count << endl;

        cout << "Outside percent = "
             << 100.0 * v.count / raw.px.size()
             << "%" << endl;

        cout << "Minimum = " << *limits.first << endl;
        cout << "Maximum = " << *limits.second << endl;
        cout << "Max violation = " << v.max_value << endl;

        Image projected = project(raw);

        Violation after = violation(projected);

        cout << "\nAFTER PROJECTION:\n";
        cout << "Outside [0,1] = " << after.count << endl;
        cout << "Max violation = " << after.max_value << endl;

        save_pgm(
            projected,
            "outputs/lab5_projected_input.pgm"
        );

        cout << "\nSaved: outputs/lab5_projected_input.pgm\n";

        return 0;
    }


    
    // LAB 5 - PROJECTED GRADIENT
    if (choice == 17) {

        cout << "\nLAB 5 - PROJECTED GRADIENT\n";

        const double lambda = calculate_lambda();

        ProjectedResult result = projected_gradient(
            data.raw,
            lambda,
            200,
            1e-6
        );

        Violation v = violation(result.image);

        cout << "Lambda = " << lambda << endl;
        cout << "Iterations = " << result.iterations << endl;
        cout << "Energy = " << result.E << endl;
        cout << "Gradient mapping = " << result.G << endl;
        cout << "Function calls = " << result.calls << endl;
        cout << "Gradient calls = " << result.grad_calls << endl;

        cout << "Total calls = "
             << result.calls + result.grad_calls
             << endl;

        cout << "Active constraints = "
             << active_constraints(result.image)
             << endl;

        cout << "Outside pixels = " << v.count << endl;
        cout << "Max violation = " << v.max_value << endl;

        cout << "MSE = "
             << mse(result.image, clean) << endl;

        cout << "PSNR = "
             << psnr(result.image, clean) << endl;

        cout << "Converged = "
             << (result.converged ? "YES" : "NO")
             << endl;

        save_pgm(
            result.image,
            "outputs/lab5_projected_gradient.pgm"
        );

        ofstream csv(
            "outputs/lab5_projected_gradient.csv"
        );

        csv << fixed << setprecision(12);

        csv << "k;E;G;active;violation;"
               "calls;grad_calls\n";

        for (const auto& row : result.rows) {

            csv << row.k << ";"
                << row.E << ";"
                << row.G << ";"
                << row.active << ";"
                << row.violation << ";"
                << row.calls << ";"
                << row.grad_calls << "\n";
        }

        csv.close();

        cout << "\nSaved: outputs/lab5_projected_gradient.csv\n";
        cout << "Saved: outputs/lab5_projected_gradient.pgm\n";

        return 0;
    }


    // LAB 5 - PENALTY METHOD
    if (choice == 18) {

        cout << "\nLAB 5 - PENALTY METHOD\n";

        const double lambda = calculate_lambda();

        PenaltyResult result = penalty_method(
            data.raw,
            lambda,
            10.0,
            7,
            100,
            1e-8
        );

        Violation v = violation(result.image);

        cout << "\nFINAL RESULTS:\n";

        cout << "Lambda = " << lambda << endl;
        cout << "Outer iterations = "
             << result.outer_iterations << endl;

        cout << "Energy = " << result.E << endl;
        cout << "Penalty P = " << result.P << endl;

        cout << "Function calls = "
             << result.calls << endl;

        cout << "Gradient calls = "
             << result.grad_calls << endl;

        cout << "Total calls = "
             << result.calls + result.grad_calls
             << endl;

        cout << "Outside pixels = "
             << v.count << endl;

        cout << "Max violation = "
             << v.max_value << endl;

        cout << "MSE = "
             << mse(result.image, clean) << endl;

        cout << "PSNR = "
             << psnr(result.image, clean) << endl;

        ofstream csv("outputs/lab5_penalty.csv");

        csv << fixed << setprecision(12);

        csv << "outer;r;Lr;h;E;P;rP;"
               "max_violation;outside;"
               "inner_iterations;calls;grad_calls\n";

        cout << "\nOUTER ITERATIONS:\n";

        for (const auto& row : result.rows) {

            cout << "s = " << row.outer
                 << " r = " << row.r
                 << " E = " << row.E
                 << " rP = " << row.rP
                 << " violation = "
                 << row.max_violation
                 << " outside = " << row.outside
                 << endl;

            csv << row.outer << ";"
                << row.r << ";"
                << row.Lr << ";"
                << row.h << ";"
                << row.E << ";"
                << row.P << ";"
                << row.rP << ";"
                << row.max_violation << ";"
                << row.outside << ";"
                << row.inner_iterations << ";"
                << row.calls << ";"
                << row.grad_calls << "\n";
        }

        csv.close();

        save_pgm(
            result.image,
            "outputs/lab5_penalty.pgm"
        );

        cout << "\nSaved: outputs/lab5_penalty.csv\n";
        cout << "Saved: outputs/lab5_penalty.pgm\n";

        return 0;
    }



if (choice == 19) {

    cout << "\nLAB 5 - EXPERIMENTS\n";

    const double lambda = calculate_lambda();
    const double L = 1.0 + 8.0 * lambda;
    const double h = 1.0 / L;

    // EXPERIMENT 1: неправильний критерій
    cout << "\nEXPERIMENT 1 - STOPPING CRITERION\n";

    ProjectedResult correct = projected_gradient(
        data.raw, lambda, 200, 1e-6
    );

    Image I = project(data.raw);
    int wrong_iterations = 0;
    double last_norm = 0.0;
    bool wrong_converged = false;

    for (int k = 0; k <= 200; k++) {

        Image g = energy_grad(I, data.raw, lambda);
        last_norm = norm(g);

        if (last_norm <= 1e-6) {
            wrong_converged = true;
            break;
        }

        if (k == 200) break;

        for (size_t i = 0; i < I.px.size(); i++) {
            I.px[i] = clamp(
                I.px[i] - h * g.px[i], 0.0, 1.0
            );
        }

        wrong_iterations++;
    }

    cout << "Correct iterations = "
         << correct.iterations << endl;

    cout << "Wrong iterations = "
         << wrong_iterations << endl;

    cout << "Wrong gradient norm = "
         << last_norm << endl;

    cout << "Wrong converged = "
         << (wrong_converged ? "YES" : "NO")
         << endl;

    // EXPERIMENT 2: множник зростання
    cout << "\nEXPERIMENT 2 - PENALTY GROWTH\n";

    ofstream growth_csv("outputs/lab5_growth.csv");
    growth_csv << fixed << setprecision(12);
    growth_csv << "growth;outer;grad_calls;"
                  "energy;violation;outside;rP\n";

    for (double growth : {2.0, 10.0, 100.0}) {

        PenaltyResult p = penalty_method(
            data.raw, lambda, growth,
            growth == 2.0 ? 40 : 7, 100, 1e-8
        );

        Violation v = violation(p.image);
        double last_rP = p.rows.back().rP;

        cout << "Growth = " << growth
             << " Outer = " << p.outer_iterations
             << " Gradient calls = " << p.grad_calls
             << " Violation = " << v.max_value
             << " rP = " << last_rP
             << endl;

        growth_csv << growth << ";"
                   << p.outer_iterations << ";"
                   << p.grad_calls << ";"
                   << p.E << ";"
                   << v.max_value << ";"
                   << v.count << ";"
                   << last_rP << "\n";
    }

    growth_csv.close();
    // EXPERIMENT 3: Lr and step versus penalty parameter.
    ofstream conditioning("outputs/lab5_conditioning.csv");
    conditioning << "outer;r;Lr;step;inner_iterations;grad_calls\n" << fixed << setprecision(12);
    PenaltyResult reference_growth = penalty_method(data.raw, lambda, 10.0, 7, 100, 1e-8);
    for (const auto& row: reference_growth.rows) {
        conditioning << row.outer << ";" << row.r << ";" << row.Lr << ";" << row.h << ";" << row.inner_iterations << ";" << row.grad_calls << "\n";
    }
    conditioning.close();

    // EXPERIMENT 4: карта активних обмежень
    cout << "\nEXPERIMENT 4 - ACTIVE CONSTRAINTS\n";

    Image active_map(
        correct.image.h, correct.image.w
    );

    int lower_count = 0;
    int upper_count = 0;
    int free_count = 0;

    for (size_t i = 0; i < correct.image.px.size(); i++) {

        double v = correct.image.px[i];

        if (v <= 1e-12) {
            active_map.px[i] = 0.0;
            lower_count++;
        }
        else if (v >= 1.0 - 1e-12) {
            active_map.px[i] = 1.0;
            upper_count++;
        }
        else {
            active_map.px[i] = 0.5;
            free_count++;
        }
    }

    save_pgm(
        active_map,
        "outputs/lab5_active_constraints.pgm"
    );

    cout << "Lower active = " << lower_count << endl;
    cout << "Upper active = " << upper_count << endl;
    cout << "Free pixels = " << free_count << endl;

    cout << "\nSaved: outputs/lab5_growth.csv\n";
    cout << "Saved: outputs/lab5_active_constraints.pgm\n";

    return 0;
}



    if (choice == 20) {
        const double lambda = calculate_lambda();
        ProjectedResult r = fista_projected(data.raw, lambda);
        Violation v = violation(r.image);
        cout << "\nLAB 5 - FISTA\nLambda = " << lambda
             << "\nIterations = " << r.iterations
             << "\nFunction calls = " << r.calls
             << "\nGradient calls = " << r.grad_calls
             << "\nTotal calls = " << r.calls + r.grad_calls
             << "\nEnergy = " << r.E
             << "\nGradient mapping = " << r.G
             << "\nPSNR = " << psnr(r.image, clean)
             << "\nOutside = " << v.count
             << "\nConverged = " << (r.converged ? "YES" : "NO") << "\n";
        save_pgm(r.image,"outputs/lab5_fista.pgm");
        ofstream file("outputs/lab5_fista.csv");
        file << "k;E;G;active;violation;calls;grad_calls\n" << fixed << setprecision(12);
        for(const auto& row : r.rows) file << row.k << ";" << row.E << ";" << row.G << ";" << row.active << ";" << row.violation << ";" << row.calls << ";" << row.grad_calls << "\n";
        return 0;
    }
    if (choice == 21) {
        cout << "\nLAB 1 - SVENN TABLE + EXPERIMENTS\n";
        ofstream tab("outputs/lab1_svenn_table.csv");
        tab << "test;k;x;f;delta;p;comment\n" << fixed << setprecision(12);
        ofstream experiments("outputs/lab1_experiments_updated.csv");
        experiments << "test;x0;delta;a;b;iterations;calls;direction\n" << fixed << setprecision(12);
        for (const auto &cfg : vector<pair<string,pair<double,double>>>{
            {"main",{X0,DELTA0}}, {"small_step",{X0,DELTA0/10}},
            {"large_step",{X0,DELTA0*10}}, {"other_start",{3.0,DELTA0}}
        }) {
            Oracle o(f); auto z = svenn(o,cfg.second.first,cfg.second.second);
            cout << cfg.first << ": [" << z.a << ", " << z.b << "] iterations=" << z.iterations << " calls=" << z.n_f << " p=" << z.p << "\n";
            experiments << cfg.first << ";" << cfg.second.first << ";" << cfg.second.second << ";" << z.a << ";" << z.b << ";" << z.iterations << ";" << z.n_f << ";" << z.p << "\n";
            for (const auto& row: z.rows) tab << cfg.first << ";" << row.k << ";" << row.x << ";" << row.f << ";" << row.delta << ";" << row.p << ";" << row.comment << "\n";
        }
        return 0;
    }
    if (choice == 22) {
        cout << "\nLAB 2 - ITERATIONS AND EPS EXPERIMENTS\n";
        ofstream tab("outputs/lab2_iterations.csv");
        tab << "eps;method;k;a;b;x1;x2;f1;f2;length\n" << fixed << setprecision(12);
        SvennResult z = svenn(oracle, X0, DELTA0);
        for (double eps: {EPS,1e-6}) for (int which=0;which<2;which++) {
            Oracle o(f);
            SearchResult res = which ? golden_section(o,z.a,z.b,eps) : dichotomy(o,z.a,z.b,eps);
            cout << "eps=" << eps << " method=" << (which?"golden":"dichotomy") << " sigma=" << res.x << " PSNR=" << 10*log10(1/res.f) << " iter=" << res.iterations << " calls=" << res.n_f << "\n";
            for (const auto& row:res.rows) tab << eps << ";" << (which?"golden":"dichotomy") << ";" << row.k << ";" << row.a << ";" << row.b << ";" << row.x1 << ";" << row.x2 << ";" << row.f1 << ";" << row.f2 << ";" << row.length << "\n";
        }
        // Demonstrate the invalid choice delta = eps by a bounded actual run.
        double a = z.a, b = z.b, eps = EPS, delta = eps;
        int stuck = 0;
        Oracle delta_oracle(f);
        for (int k=0;k<40 && b-a>eps;++k) {
            double x1=(a+b-delta)/2.0, x2=(a+b+delta)/2.0;
            double y1=delta_oracle.value(x1), y2=delta_oracle.value(x2);
            if (y1<=y2) b=x2; else a=x1;
            ++stuck;
        }
        cout << "Experiment delta=eps: iterations=" << stuck << " length=" << b-a
             << " met_tolerance=" << (b-a<=eps ? "YES":"NO") << "\n";
        return 0;
    }
    if (choice == 23) {
        const double lambda = calculate_lambda();
        Image sol = tikhonov_solve(data.raw,lambda);
        Violation v = violation(sol);
        cout << "\nLAB 5 - UNCONSTRAINED RAW BASELINE\nLambda=" << lambda << "\nOutside=" << v.count << "\nMax violation=" << v.max_value << "\nEnergy=" << energy(sol,data.raw,lambda) << "\n";
        save_pgm(sol,"outputs/lab5_unconstrained.pgm");
        return 0;
    }

    cout << "Unknown option.\n";

    return 0;
}
