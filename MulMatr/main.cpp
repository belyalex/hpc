#include <random>

#include "matrix.h"
#include "timer.h"

int main() {
    std::default_random_engine generator(time(nullptr));
    std::uniform_real_distribution<double> distribution;

    std::cout << "Matrix size\tShtrassen mul\tSimple mul\tEquals\n";

    for (int n = 512; n <= 4096; n *= 2) {
        const int D = (n < 4096 ? 5 : 1);
        matrix a(n), b(n);

        for (size_t i = 0; i < a.size(); i++) {
            for (size_t j = 0; j < a.size(); j++) {
                a[i][j] = distribution(generator);
                b[i][j] = distribution(generator);
            }
        }

        std::cout << n << '\t';

        matrix c;
        {
            timer t(D);
            for (int i = 0; i < D; i++) {
                c = shtrassen_mul(a, b);
            }
        }
        std::cout << '\t';
        matrix d;
        {
            timer t(D);
            for (int i = 0; i < D; i++) {
                d = a * b;
            }

        }
        std::cout << '\t';
        std::cout << ((c == d) ? "true" : "false") << '\n';
    }
    return 0;
}
