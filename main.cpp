#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"
#include "display_tools.h"

int main()
{
    vector<Matrix> localKs;
    vector<Vector> localbs;

    Matrix K;
    Vector b;
    Vector T;

    cout << "IMPLEMENTACION DEL METODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- ECUACIONES DE NAVIER-STOKES\n"
         << "\t- 1 DIMENSION\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n"
         << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m);

    crearSistemasLocales(m, localKs, localbs);

    zeroes(K, m.getSize(NODES) * 2);
    zeroes(b, m.getSize(NODES) * 2);

    ensamblaje(m, localKs, localbs, K, b);

    applyDirichlet(m, K, b);

    showMatrix(K);

    // zeroes(T, b.size());

    // calculate(K, b, T);

    // cout << "LA RESPUESTA ES: \n";

    // showVector(T);

    return 0;
}
