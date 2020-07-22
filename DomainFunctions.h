#ifndef GUARD_DOMAINFUNCTIONS_H
#define GUARD_DOMAINFUNCTIONS_H
#include "main.h"
void UpdateDomain(cell**&, void function(cell&));
void UpdateGhostCells(cell**&, int&);
void UpdateMainCells(cell**&, int&, void function(cell&, int&));
void UpdateMainCells(cell**&, int&, double&, void function(cell&, int&, double&));
void UpdateMainCells(cell**&, void function(cell&));
void CalculateFlux(cell**, int&);
void FirstOrder(cell&, int&);
void WENO(cell&, int&);
void RungeKutta(cell&, int&, double& dt);
void RungeKuttaReset(cell&);
void ZeroFlux(cell&, int&);
void Initialize(cell**&);
void InitialCondition(cell**&);
void Delete(cell&);
double& MinRequiredTimeStep(cell**&, double&);
double CellTimeStep(cell&);
double* F(double*&);
double* G(double*&);
#endif /* DOMAINFUNCTIONS_H */

