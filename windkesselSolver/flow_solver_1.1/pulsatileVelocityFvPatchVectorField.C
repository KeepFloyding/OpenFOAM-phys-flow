/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "pulsatileVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "pulsatileVelocityFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pulsatileVelocityFvPatchVectorField::pulsatileVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0)
{}


pulsatileVelocityFvPatchVectorField::pulsatileVelocityFvPatchVectorField
(
    const pulsatileVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


pulsatileVelocityFvPatchVectorField::pulsatileVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("pulsatileVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    evaluate();
}


pulsatileVelocityFvPatchVectorField::pulsatileVelocityFvPatchVectorField
(
    const pulsatileVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    maxValue_(fcvpvf.maxValue_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pulsatileVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * *  Pulsatile Calculation  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Possible culprits
// Wrong R value

   
//------------------------------------------------------------------------------   
// Define all constants and parameters
//------------------------------------------------------------------------------   

	scalar pi = constant::mathematical::pi;

	boundBox bb(patch().patch().localPoints(), true);

	scalar R = ((bb.max() - bb.min()) & y_)/2;
	scalar mu  = 3.5E-3; // Dynamic Viscosity (Pa s)
	scalar T   = 1.0;    // Period (s) 
	scalar rho = 1060;   // Blood Density (kg/m^3)
	scalar omega = 2*pi/T;            	// Angular Velocity (1/s)
	scalar alpha = R*pow((omega*rho/mu),0.5);	// Womersely Number (dimensionless)

	// Middle co-ordinates of inlet 

    // Get range and orientation

	const scalar time = 0.1007;    			// Current time (s)
	scalar Q;
	
	// For the Bessel functions 
	scalar a_1, b_1, Ber_prev, Bei_prev, B_1r_prev, B_1i_prev,  fact_2n, fact_2n_p_1;
	scalar Ber, Bei, B_1r, B_1i;
	scalar a_3, b_3, a_4, b_4, a_5, b_5, denom_r, denom_i;
	scalarField nom_r, nom_i, a_2, b_2;
	
	// For the loops
	scalar n, n_2, it, it_2;

//------------------------------------------------------------------------------   
// Finding the flowrate at the given time; coefficients and expression 
//------------------------------------------------------------------------------   
	
	scalar a0 =      4.307  ;
	scalar a1 =     -0.2427, b1 =      1.125   ;
	scalar a2 =     -0.9454, b2 =      0.5766  ;
	scalar a3 =     -0.4395, b3 =     -0.3121  ;
	scalar a4 =     -0.2677, b4 =    -0.4183  ;
	scalar a5 =      0.3666, b5 =    -0.2777  ;
	scalar a6 =     0.07413, b6 =    0.05799  ;
	scalar a7 =     0.111,   b7 =      -0.01202  ;
	scalar a8 =     0.0431,  b8 =     0.1057  ;
	scalar w  =       2*pi  ;
	

	
	Q = (a0 + a1*cos(time*w) + b1*sin(time*w) + a2*cos(2*time*w) + b2*sin(2*time*w) + a3*cos(3*time*w) + b3*sin(3*time*w) + 
		a4*cos(4*time*w) + b4*sin(4*time*w) + a5*cos(5*time*w) + b5*sin(5*time*w) + a6*cos(6*time*w) + b6*sin(6*time*w) + a7*cos(7*time*w) 
		+ b7*sin(7*time*w) + a8*cos(8*time*w) + b8*sin(8*time*w))*1E-6; //UNITS ARE IN M^3/S
		
		
//------------------------------------------------------------------------------   	
// Evaluating the constant Bessel functions, J_0(alpha*i^(3/2)) and J_1(alpha*i^(3/2))
//------------------------------------------------------------------------------   

	// Real and complex expression for i^(3/2)
	a_1 = cos(3*pi/4);
	b_1 = sin(3*pi/4);
	
	// Real and complex components of J_0 
	Ber_prev = 0;
	Bei_prev = 0;
	
	// Real and complex components of J_1 
	B_1r_prev = 0;
	B_1i_prev = 0;
	

//------------------------------------------------------------------------------   
// Calculation of factorials and Bessel constants
//------------------------------------------------------------------------------   
	 
	
	fact_2n_p_1 = 1; // (2n + 1)!
	
	for(n = 0; n < 11;n++)
	{
		fact_2n = 1;    // 2n!
		if (n>0)
        {
			for (it = 1; it <= 2*n;it++)
			{
				fact_2n = fact_2n*it;
			}
        

        }    
   		fact_2n_p_1 =  fact_2n*(2*n+1);


    Ber = Ber_prev + pow(-1,3*n)*pow((alpha/2),4*n)/(fact_2n*fact_2n);
    Ber_prev = Ber;
    
    Bei = Bei_prev + pow(alpha/2,(4*n+2))*pow(-1,3*n)/(fact_2n_p_1*fact_2n_p_1);
    Bei_prev = Bei;

    B_1r = B_1r_prev + a_1*pow((alpha/2),(4*n+1))*pow(-1,3*n)/((fact_2n*fact_2n)*(2*n +1)) - b_1*pow((alpha/2),(4*n+3))*pow(-1,3*n)/((fact_2n_p_1*fact_2n_p_1)*(2*n +2));
    B_1r_prev = B_1r;
    
    B_1i = B_1i_prev + b_1*pow((alpha/2),(4*n+1))*pow(-1,3*n)/((fact_2n*fact_2n)*(2*n +1)) + a_1*pow((alpha/2),(4*n+3))*pow(-1,3*n)/((fact_2n_p_1*fact_2n_p_1)*(2*n +2));
    B_1i_prev = B_1i;
    
	}
	
//------------------------------------------------------------------------------   	 
// Finding the value of r 
//------------------------------------------------------------------------------   
	 
    
    // Finding the centre of the domain
    
   
	vector ctr = 0.5*(bb.max() + bb.min());
       
    // Evaluating the centre of each face
       
    const vectorField& c = patch().Cf();
    
    // Determing the value of the distance r
	scalarField r = pow((c - ctr) & (c - ctr),0.5);


	forAll(r,iter)
	{
		
		if (r[iter] >R)
		{
			r[iter] = R;

		}
	}
    
//------------------------------------------------------------------------------   	     
// Calculation of the value of the Bessel function for different r
//------------------------------------------------------------------------------   
	        
 		scalarField Ber1_prev = ((c - ctr) & y_)*0;
        scalarField Bei1_prev = ((c - ctr) & y_)*0;
        
        scalarField Ber_1, Bei_1;
        
		for (n_2 = 0; n_2<10; n_2++)
		{
            // Calculating the factorials 
			fact_2n = 1;
            if (n_2 > 0)
            {   
                for (it_2 = 1; it_2 <= 2*n_2; it_2++)
				{
                    fact_2n = fact_2n*it_2;
                }
            }
            
            fact_2n_p_1 =  fact_2n*(2*n_2+1);


			Ber_1 = Ber1_prev + pow(-1,3*n_2)*pow((alpha*r/(2*R)),4*n_2)/(fact_2n*fact_2n);
			
			
			Bei_1 = Bei1_prev + pow((alpha*r/(2*R)),(4*n_2+2))*pow(-1,3*n_2)/(fact_2n_p_1*fact_2n_p_1);
        
			Ber1_prev = Ber_1;
			Bei1_prev = Bei_1;


        }    
		
//------------------------------------------------------------------------------   	
// Calculation of the velocity for r   
//------------------------------------------------------------------------------   
		
		a_2 = Ber_1;
        b_2 = Bei_1;
        
        nom_r = 1 - (Ber*a_2 + b_2*Bei)/(Ber*Ber + Bei*Bei);
        nom_i = (a_2*Bei - b_2*Ber)/(Ber*Ber + Bei*Bei);
        
        a_3 = B_1r; 
		b_3 = B_1i;
        a_4 = Ber; 
		b_4 = Bei;
        a_5 = a_4*a_1 - b_4*b_1;
		b_5 = b_4*a_1 + b_1*a_4;
        
        denom_r = 1 - (2/alpha)*(a_3*a_5+b_5*b_3)/(a_5*a_5 + b_5*b_5); 
        denom_i = -(2/alpha)*(b_3*a_5 - b_5*a_3)/(a_5*a_5 + b_5*b_5);
        
//------------------------------------------------------------------------------   	
// Assigning the velocity to the BC  
//------------------------------------------------------------------------------   
	
	
	
    // Pulsatile Velocity Profile	

	vectorField vel_r = n_*(Q/(pi*R*R))*((nom_r*denom_r + nom_i*denom_i)/(denom_r*denom_r + denom_i*denom_i));

    	vectorField::operator = (vel_r);
    

    
}


// Write
void pulsatileVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue")
        << maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, pulsatileVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
