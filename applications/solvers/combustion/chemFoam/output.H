    runTime.write();

    Info<< "Qdot = " << Qdot
        << ", T = " << thermo.T()[0]
        << ", p = " << thermo.p()[0]
        << ", " << Y[0].name() << " = " << Y[0][0]
        << endl;

    post<< runTime.value() << token::TAB << thermo.T()[0] << token::TAB
        << thermo.p()[0] << token::TAB << Qdot << token::TAB << integratedHeat << token::TAB << rho[0];

	forAll(Y,i)
        {
            post<< token::TAB << Y[i][0];
        }
	post<< endl;
