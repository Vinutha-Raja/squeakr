//
// Created by Vinutha Raja on 3/28/23.
//

#include "minimizer.h"
#include <vector>
using namespace std;

int main() {
    string seq = "ACGATCGACG";
    MinimizerScanner scanner(10, 5, 0);
    scanner.LoadSequence(seq);
    uint64_t *mmp;
    while ((mmp = scanner.NextMinimizer()) != nullptr)
        printf("%016lx\n", *mmp);
    return 0;
}
