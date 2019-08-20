dtChem = chemistry.solve(runTime.deltaT().value());
scalar Qdot = chemistry.Qdot()()[0]/rho[0];
integratedHeat += Qdot*runTime.deltaT().value();
//totalIntHeat += chemistry.Qdot()()[0]*runTime.deltaT().value();
