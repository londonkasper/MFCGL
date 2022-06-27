// MFCGLView.cpp : implementation of the CMFCGLView class
//

#include "afx.h"
#include "stdafx.h"
#include "MFCGL.h"

#include "MFCGLDoc.h"
#include "MFCGLView.h"
#include <omp.h> 
#include <ctime>

#include "Math.h"
#include <stdlib.h>
#include <stdio.h>
//#include <iostream>
#include <fstream>

//using namespace std;


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMFCGLView

IMPLEMENT_DYNCREATE(CMFCGLView, CView)

BEGIN_MESSAGE_MAP(CMFCGLView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
	ON_WM_SIZE()
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_RBUTTONDOWN()
	ON_WM_TIMER()
	ON_COMMAND(ID_AUTO_RUN, &CMFCGLView::OnAutoRun)
	ON_COMMAND(ID_AUTO_STOP, &CMFCGLView::OnAutoStop)
	ON_COMMAND(ID_AUTO_SINGLECPU, &CMFCGLView::OnAutoSinglecpu)
	ON_COMMAND(ID_AUTO_ADDCPU, &CMFCGLView::OnAutoAddcpu)
	ON_COMMAND(ID_AUTO_MAXCPU, &CMFCGLView::OnAutoMaxcpu)
	ON_COMMAND(ID_AUTO_REMOVECPU, &CMFCGLView::OnAutoRemovecpu)
	ON_WM_SHOWWINDOW()
	ON_COMMAND(ID_MOVE, &CMFCGLView::OnMove)
	ON_COMMAND(ID_ShowMolecule, &CMFCGLView::OnShowmolecule)
	ON_COMMAND(ID_NEWSIM, &CMFCGLView::OnNewsim)
	ON_COMMAND(ID_WIRE, &CMFCGLView::OnWire)
	ON_COMMAND(ID_RUN, &CMFCGLView::OnRun)
	ON_COMMAND(SHEAR, &CMFCGLView::OnShear)
	ON_COMMAND(PERIODIC, &CMFCGLView::OnPeriodic)
	ON_COMMAND(ID_PRESSURE, &CMFCGLView::OnPressure)
	ON_COMMAND(ID_INIT, &CMFCGLView::OnInit)
	ON_COMMAND(ID_PROPDLG, &CMFCGLView::OnPropdlg)
	ON_COMMAND(ID_AUTO_TIMER, &CMFCGLView::OnAutoTimer)
	ON_COMMAND(ID_DAMP1, &CMFCGLView::OnDamp1)
	ON_COMMAND(ID_DAMPDN, &CMFCGLView::OnDampdn)
	ON_COMMAND(ID_DAMPUP, &CMFCGLView::OnDampup)
	ON_COMMAND(ID_FILE_OPEN, &CMFCGLView::OnFileOpen)
	ON_COMMAND(ID_FILE_SAVE, &CMFCGLView::OnFileSave)
	ON_COMMAND(ID_WallDown, &CMFCGLView::OnWalldown)
	ON_COMMAND(ID_REVERSEFLOT, &CMFCGLView::OnReverseflot)
	ON_COMMAND(ID_BOLTZMANN_VEL, &CMFCGLView::OnBoltzmannVel)
	ON_COMMAND(ID_OPENWALL, &CMFCGLView::OnOpenwall)
	ON_COMMAND(ID_MSAVE, &CMFCGLView::OnMsave)
	ON_COMMAND(ID_Post, &CMFCGLView::OnPost)
	ON_WM_MOUSEWHEEL()
	ON_COMMAND(ID_VISUALCLOSE, &CMFCGLView::visualclose)
	ON_COMMAND(ID_LATTICE, &CMFCGLView::OnLattice)
	ON_COMMAND(ID_S, &CMFCGLView::OnS)
END_MESSAGE_MAP()
  
// CMFCGLView construction/destruction
 
CMFCGLView::CMFCGLView()
{
////MB////
visual = 0;
Nfolder = 1;
////MB////

	m_strTemp = "";
	m_strNM   ="";
	m_Gcount = 0;
	m_rf = 1.0;

	///////////////////////////////////////
	// Program    Property               //
	///////////////////////////////////////

	//5    2.5    9.0     5nm  ANGLE  50
	//-8.5  -5.0    -8.0   10nm  ANGLE  80
	//-18.5 -10.0 -8.0  20nm ANGLE 120
	//-18.5 -25.0 -8.0  50nm ANGLE 147
	view_angle = 50.0;
	m_xRotate=0.0;
	m_yRotate=0.0;
	m_zRotate=0.0;
	//Set the origin of the coordinate 
	m_xTranslate =-8.55f;
	m_yTranslate =-4.5f;
	m_zTranslate =15.0; //왜 이게 6.0코드와는 반대로 되어있을까?  m_zTranslate =-15.0;
	m_dautocolor = 0.5;

	//Flags 1: ON  2:OFF 
	m_nMode = 1;           //For rotation of the cooridnate
	m_nShow_molecule = 0;  //[for draw 1,2,3,4th moleculr KE and PE
	m_nWire = 2;           //Show 
	m_nPeriodic = 1;       //Show Image for periodic bc
	m_nSp = 0;             //0:off, 1: Pressure driven case, 2:Shear driven case

	m_ncutoff = 3;         //Cutoff for Smart wall wall
	m_nWrange = m_ncutoff; //For smart wall range
	m_nwallinteraction = 1;//Including wall intraction
	m_time = 0;            //Timer flag (use animation)
	m_n3d = 0;             // 3D flag [NOT USE]
	m_nroop = 5;           // number of loop per screen update;
	m_nmilisec = 1;       // Call timer interval
	m_nWallDown = 0;       //Move Wall Down

	d_ntimestep = 0.0;    //for recording timestep

	m_istrt = 1.0f;
	m_acc = 0.1f;
	m_max = 0.5f;
	m_spring = 3.0f;
	m_fH2 = 6.0f;  
	m_fW2 = 8.0f;  
	m_fZ2 = 4.0f;  

	///////////////////////////////////////
	// Simulation Property               //
	///////////////////////////////////////

	//Reference/////////////////////////////////////////////////////////
	
	//[J = kg m2 s-2]
	//tera = e12	giga = e9	mega =e6	kilo=e3	   hecto=e2  deca=e1 	
	//deci= e-1		centi=e-2	milli=e-3	micro=e-6  nano=e-9	 pico=e-12
	//femto=e-15	atto=e-18	zepto=e-21  yocto=e-24

	//Unit conversion
	//first unit * this = second unit, e.g.) 1m * m_to_nm = 1e9nm
	//second unit/this = first unit
	m_to_nm = 1e9;
	s_to_ps = 1e12;
	kg_to_g = 1e3;

	joule_to_gnmps = (kg_to_g*m_to_nm*m_to_nm)/(s_to_ps*s_to_ps); //= 1e-3;
	//multiple when change [J = kg m2 s-2] to [g nm2 ps-2]
	 //divide when change [g nm2 ps-2] to [J = kg m2 s-2]

	//Reference/////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////
	//All Units used in the simulation is [g],[nm],[ps],[g/mole],[K]//
	// All constant defined here and MUST convert to gnmpsK unit	//
	//////////////////////////////////////////////////////////////////

	m_AvogN = 6.02214199e23;   //Avogadro's number      [atoms] per [mole]  Non-dimensional number
	m_UGasN = 8.3144;         //Universal Gas Constant [J][mol-1][K-1]
	m_UGasN = 8.3144*joule_to_gnmps; //to gnmps unit

	m_BoltzN = 1.3806503E-23;     //Boltzmann's Constant   [J][K-1]
	m_BoltzN = 1.3806503E-23*joule_to_gnmps;

	//check if the number is correct
	m_UGasN = m_BoltzN*m_AvogN;  //Universal Gas Constant [J][mol-1][K-1]

	//[these value is reset at below]	m_BoltzN = 1.3806503E-23;     //Boltzmann's Constant   [J][K-1]
	//for argon.  // Not used in the Code. because the epsylon is already divided by mass

	//Based on [argon]
	m_dMass_wall = 39.948;     //Atomic Mass [g/mole] - wall molecules
	m_dMass_gas = 39.948;      //Atomic Mass [g/mole]- gas molecules 

	m_massofargon = m_dMass_gas/m_AvogN; // [gram] = 6.6335201106741e-023g 
	m_Kb = m_massofargon/m_BoltzN;     // m/Kb  Mass divided by Boltzmann constant in nmps unit [argon]       
	m_massofargon = m_Kb*m_BoltzN; //6.6335201106741e-023[g]
	
	m_dEgas = 119.8*m_BoltzN/m_massofargon; //= 0.024934258610877

	double temp;

	temp = 120*1.3806503E-23*m_AvogN;
	temp = temp*2.390e-4;
	
	m_dEgas = 4.0*m_dEgas;  //Epsylon - gas molecules   [nm]2/[ps]2 //Set this to e/m
	m_dEwall = m_dEgas*1.0; //Epsylon - wall molecules  //Set this to e/m 
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//m_fMdist = 0.381868235511538;       //[nm]//Distance between Wall Molecules [nm]
	//m_fMdist = 0.382; // equilibrium distance
	m_fMdist = 0.540; // for lattice parameter  
	m_dSigma = 0.3405;        //[nm]//Sigma(Diameter) in LJ poetntial eqn
	m_fMrad = m_dSigma/2.0; //[nm]//Radius of Sphere [nm]

	//m_dGravP = 0.03;        //Gravity term for flow drive
	double tau,ChaT;
	tau = m_dSigma*sqrt(4.0/m_dEgas); //= 2.1563481794417
	ChaT = (m_dEgas/4.0)*m_Kb; //e/m * m/kb  = 119.8
	m_dGravP = 0.1*m_dSigma/(tau*tau);        //Gravity term for flow drive 	//m_dGravP = 0.0072433706696153
	
///// 1.6539588000000e-021 [J] for epsylon  kb * 119.8K;
	//m_Umax = (4.0)*(0.5)*0.5*sqrt(m_dEgas/4.0); //for 4h same shear rate need 0.5 for shear both direction
	//m_Umax = (2.0)*(0.5)*0.5*sqrt(m_dEgas/4.0); //for 2h same shear rate
	m_Umax = 0.5*sqrt(m_dEgas/4.0); //for 1h same shear rate
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//m_Umax = (0.75)*0.5*sqrt(m_dEgas/4.0); //for 0.75U for 3rd work different shear rate.
	//Pradip Velocity
	//m_Umax = 0.05*sqrt((5.0/3.0)*m_Ttarget/m_Kb);  //fluid sumulation that I use for reproduce pradip dsmc
	//m_Umax = 0.5*sqrt(m_dEgas/4.0);  //fluid simulation that I use
	//m_aTwM[ii].m_fX =m_aTwM[ii].m_fX + m_Umax*m_fdt;// m_dmove*3;
	//m_Umax = sqrt(m_dEgas/4.0);
    //m_Umax = 0.15704673549628 nm/ps
	//m_Umax = 0.5*sqrt(m_dEgas/4.0);
	//m_Umax = 0.078523367748142

	m_nwMs = 30000;         //Number of Wall Molecules
	m_nMs = 6300;           //Number of Inner Molecules
	m_fW = 0.540*m_fW2*1.0;	    //[nm]   //Channel Length
	m_fH = 0.540*m_fH2*1.0;       //[nm]   //Width of the Channel 
	m_fZ = 0.540*m_fZ2*1.0;       //[nm]   //depth of the Channel 


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//	4.32000000000
//	3.24000000000
//	2.16000000000 

	//Number Density
	d_numden = (578.0/(m_fW*(m_fH-m_dSigma*1.0)*m_fZ)) * (m_dSigma*m_dSigma*m_dSigma);
	//Mass density
	d_numden = (578.0*m_massofargon/(m_fW*(m_fH-m_dSigma*1.0)*m_fZ))*1e21 ;
	//1.2682047642535 [g cm-3] g cm unit for ref comparison.

	m_fdt = 0.004;           //[ps]//Time step - 0.01ps = 10fs 
	m_Ttarget = 100.0;        //[K] Target temperature 
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	m_ndamping = 0;        //Thermostat 0:Normal user controll, 1:velocity scalling 2:Nose hoover 
	m_dampcoef = 1;        //scale velocity

	m_gamma = 0.0;  //for shear 
	m_gamma2 = 0.0;  // for pressure driven
	//m_gamma = 2.0*m_Umax/m_fH; defined on OnShear()

	//1.783 kg/m^3 argon density STP=298K 1atm 
	//64.41m mean free path.

	// sqrt(pi)/2*Kn -> 1.7724538509055158/2*64.41/50
	// -> 1.14163752536824272678 k=1 case
	//1.783 / m_massofargon = 26639165759287257496175580.968866
	//cubic root of this -> 298657578.8644529 per 1m.
	//298.6575788644529 per 1micro meter
	//29.86575788644529 per 100nm
	//14.932878943222645 per 50nm
	//(14.932878943222645^2)*3 = 668.97262059882677674073412238808
	// -> 670/10 = 67
	//0.2986575788644529 per 1nano meter
    // about 6 molecules for this qubic.

	//Target temperature (setting) 
	wallTempT=100.0;         
	wallTempB=100.0;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//d_K_M = 5.0*16*4.0*m_dEgas/(m_dSigma*m_dSigma);
	//d_K_M = 2.0*16*4.0*m_dEgas/(m_dSigma*m_dSigma);
	//d_K_M = 2.0*16.0*4.0*m_dEgas/(m_dSigma*m_dSigma);  //2k
	d_K_M = m_spring*16.0*4.0*m_dEgas/(m_dSigma*m_dSigma);  //3k
	//Spring constant for solid wall.
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	d_vibw = 0.0;
 
	m_nWtype = 1;           //Wall type 0: BCC  1:FCC
	m_nNrow = 2;            //Number of the wall layer
	//m_dcutdist = 2.5*m_dSigma;       //[nm]//LJ cut off distance
	m_cutpot = 50000.0;     //Cut off Interaction force (Prevent sudden divergence)
	//m_cellsize = 1.0; // cell for linked_cell algorithm
	m_cellsize = 0.540*2.0;
	m_dcutdist = m_cellsize;       //[nm]//LJ cut off distance

//need to fix corresponding force switch.

	m_ygsize = 0.10;  // size for y direction 
	m_bxsize = 0.5; //
	m_bysize = 0.5; //

	///////////////////////////////////////
	// Initialization                    //
	///////////////////////////////////////

	m_aTwM = new CMolecular [m_nwMs];  //Top wall molecules
	m_aBwM = new CMolecular [m_nwMs];  //Bottom wall molecules

	m_aTwM_OMP = new CMolecular [m_nwMs];  //Temporary storage for OpenMP
	m_aBwM_OMP = new CMolecular [m_nwMs];  //Temporary storage for OpenMP


	m_aTwM_orig = new CMolecular [m_nwMs];  //Top wall molecules Original position
	m_aBwM_orig = new CMolecular [m_nwMs];  //Bottom wall molecules Original position
	m_aTwM_old = new CMolecular [m_nwMs];  //Top wall molecules Original position
	m_aBwM_old = new CMolecular [m_nwMs];  //Bottom wall molecules Original position

	m_aM_OMP = new CMolecular [m_nMs];     //Temporary storage for OpenMP
	m_aM = new CMolecular [m_nMs];     //Gas molecules


	m_V  = 0.0;
 
	//variables for system property
	m_T  = 0.0;
	m_P = 0.0;
	m_Tsc = 0.0;
	
	//for wall generator
	m_nz = 1;
	m_nx = 1;
	m_ny = 1;
	
	m_nnz = 1;
	m_nnx = 1;
	m_nny = 1;


	//For visualization
	d_Twalltemp = 0.0; //variable for Wall Temperature
	d_Bwalltemp = 0.0; //variable for Wall Temperature

	d_Twalltemp1 = 0.0; //variable for Wall Temperature
	d_Bwalltemp1 = 0.0; //variable for Wall Temperature
	d_Twalltemp2 = 0.0; //variable for Wall Temperature
	d_Bwalltemp2 = 0.0; //variable for Wall Temperature
	d_Twalltemp3 = 0.0; //variable for Wall Temperature
	d_Bwalltemp3 = 0.0; //variable for Wall Temperature
	
	m_Kin = 0.0;
	m_Pot = 0.0;
	m_Tot = 0.0;

	m_Potx = 0.0;
	m_Poty = 0.0;

	m_Kinsc = 0.0;
	m_Potsc = 0.0;
	m_Totsc = 0.0;

	//For Nose-Hoover Thermostat;
	m_zeta = 0.0;
	m_zeta_dot = 0.0;
	m_tau = 0.0;
	m_Q = 0.0;

	
	m_heatfluxx = 0.0;
	m_heatfluxy = 0.0;
	m_heatfluxz = 0.0;
	m_heatsumy = 0.0;
	m_heatsumx = 0.0;
	m_heatsumz = 0.0;



	m_dampcoef = 1.0;

	int ii;
		for(ii = 0; ii<m_nMs ; ii++)
		{
		m_aM[ii].m_fStress = 0.0; //initialize
		m_aM[ii].m_fVirialxy = 0.0; //initialize
		}

		double dNCX,dNCY,dNCZ;
		dNCX = m_fW/m_cellsize;
		dNCY = m_fH/m_cellsize;
		dNCZ = m_fZ/m_cellsize;


	m_vec.fsignx = 0.0;	m_vec.fsigny = 0.0;
	m_vec.accx = 0.0;	m_vec.accy = 0.0;	m_vec.fForce = 0.0; 
	m_vec.fPot = 0.0;


	m_nMs = 0;				//Initialize number of wall molecule
	m_nwMs = 0;            //Initialize number of wall molecule

		for (int jj = 0; jj < 10 ; jj++)
		{	
			strG[jj] = "0";	
			strGDen[jj] = "0";

			dG[jj] = 0.0;	
			denCount[jj] = .0;
		}

	//For data collection and time averaging!
	for (ii =0; ii<500 ; ii++)
	{
	 m_data_count5[ii] = 0.0;	 m_data_count10[ii] = 0.0;
	 m_data_count50[ii] = 0.0;	 m_data_count100[ii] = 0.0;
	 m_data_count500[ii] = 0.0;

	 m_data_temp5[ii] = 0.0;	 m_data_temp10[ii] = 0.0;
	 m_data_temp50[ii] = 0.0;	 m_data_temp100[ii] = 0.0;
	 m_data_temp500[ii] = 0.0;

	 m_data_collectT5[ii] = 0.0;	 m_data_collectT10[ii] = 0.0;
	 m_data_collectT50[ii] = 0.0;	 m_data_collectT100[ii] = 0.0;
	 m_data_collectT500[ii] = 0.0;

	 m_data_collectD5[ii] = 0.0;	 m_data_collectD10[ii] = 0.0;
	 m_data_collectD50[ii] = 0.0;	 m_data_collectD100[ii] = 0.0;
	 m_data_collectD500[ii] = 0.0;

	}



	 m_startstep1 = 11000*3;
	 m_endstep1 = 30000*3;
	 m_averagingstep1 = m_endstep1 - m_startstep1;

	 m_startstep2 = 31000*3;
	 m_endstep2 = 50000*3;
	 m_averagingstep2 = m_endstep2 - m_startstep2;
	
	 m_startstep3 = 51000*3;
	 m_endstep3 = 70000*3;
	 m_averagingstep3 = m_endstep3 - m_startstep3;
	
	 m_startstep4 = 71000*3;
	 m_endstep4 = 90000*3;
	 m_averagingstep4 = m_endstep4 - m_startstep4;
	
	 m_startstep5 = 91000*3;
	 m_endstep5 = 110000*3;
	 m_averagingstep5 = m_endstep5 - m_startstep5;
	//For data collection and time averaging!

	 d_totalLoop = 0.0;

	//Heat
	for (ii =0; ii<100000 ; ii++)
	{
	d_J[ii][0]=0.0;
	d_J[ii][1]=0.0;
	d_J[ii][2]=0.0;
	}
	n_Jcount = 0;
	n_Jtotal = 0;



	 
	int xx,yy,zz;
	//m_dcutdist
	////For periodic boundary condition and sort for cell

	
		for (xx = 0; xx <= int(dNCX)+1 ;xx++) //Initialize cell data
		{	
			for (yy = 0; yy <= int(dNCY)+1 ;yy++)
			{
				for (zz = 0; zz <= int(dNCZ)+1 ;zz++)
				{
					Cell[xx][yy][zz].NM = 0;
					for (ii = 0; ii < 64 ;ii++)
						{
						Cell[xx][yy][zz].N[ii] = 0;
						Cell[xx][yy][zz].x[ii] = 0.0;
						Cell[xx][yy][zz].y[ii] = 0.0;
						Cell[xx][yy][zz].z[ii] = 0.0;
						Cell[xx][yy][zz].vx[ii] = 0.0;
						Cell[xx][yy][zz].vy[ii] = 0.0;
						Cell[xx][yy][zz].vz[ii] = 0.0;
						}


				}
			}
		}




/*
주기율표 제0족에 속하는 비활성기체원소. 
 
   원소기호  Ar 
   원자번호  18 
   원자량  39.948 
   녹는점  -189.2℃ 84K
   끓는점  -185.7℃ 89K
   비중  1.7834g/ℓ 
 
 본문 

 1894년에 영국의 J.W.레일리가 W.램지의 협력을 얻어 공기에서
 산소를 제거하고 얻은 질소의 비중과, 질소화합물을 분해하여
 얻은 질소의 비중이 다른 점에 착안하여, 공기에서 이 물질을
 분리시켜 발견하였다.  이 물질은 그때까지 알려진 원소와 달리
 화학적 성질이 극히 활발하지 않아 모든 물질과 반응하지 않았던
 데서, 그리스어인 'an ergon(게으름쟁이)'을 따서 명명되었다.
 */




	// TODO: add construction code here
	nRotate =0;  //0=fix, 1=Rotate,
	iCPU = 1;
	glcrossx = 0.0f;
	glcrossy = 0.0f;

	duration = 0.0;
	start =0.0;
	finish = 0.0;


}

CMFCGLView::~CMFCGLView()
{
	KillFont();

	delete [] m_aTwM;
	delete [] m_aBwM;

	delete [] m_aTwM_OMP;
	delete [] m_aBwM_OMP;

	delete [] m_aTwM_orig;
	delete [] m_aBwM_orig;
	delete [] m_aTwM_old;
	delete [] m_aBwM_old;

	delete [] m_aM;
	delete [] m_aM_OMP;


}

BOOL CMFCGLView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CMFCGLView drawing

void CMFCGLView::OnDraw(CDC* /*pDC*/)
{
	CMFCGLDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	//wglMakeCurrent(m_hDC, m_hRC);

	DrawScene();

	SwapBuffers(m_pDC->GetSafeHdc());

	//wglMakeCurrent(m_hDC, NULL);

	// TODO: add draw code for native data here
}


// CMFCGLView printing

BOOL CMFCGLView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CMFCGLView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CMFCGLView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}


// CMFCGLView diagnostics

#ifdef _DEBUG
void CMFCGLView::AssertValid() const
{
	CView::AssertValid();
}

void CMFCGLView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMFCGLDoc* CMFCGLView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMFCGLDoc)));
	return (CMFCGLDoc*)m_pDocument;
}
#endif //_DEBUG


// CMFCGLView message handlers

void CMFCGLView::OnSize(UINT nType, int cx, int cy)
{
//	cx = 2.0*cy;
	CView::OnSize(nType, cx, cy);

	GLdouble dAspect = (GLdouble) cx/ (GLdouble) cy;
	glMatrixMode(GL_PROJECTION);  // Projection Matrix Mode
	glLoadIdentity();             // Initilaze Matrix
	gluPerspective(50.0f, dAspect, 1.0, 40.0);   // Viewing volume //BHK VOLUME
	glViewport(0, 0, cx, cy);     // Set the size of viewport to entire domain

	// TODO: Add your message handler code here

	//VERIFY(wglMakeCurrent(m_hDC, m_hRC));
	
	//GLResize(cx,cy);

	//VERIFY(wglMakeCurrent(NULL, NULL));
	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());


}

void CMFCGLView::GLResize(int cx, int cy)
{

	GLfloat fAspect;

	if (cy == 0)
		cy = 1;

	glViewport(0,0,cx,cy); //, 모니터 전체 화면을 뷰포트 영역으로 잡아 준다

	fAspect = (GLfloat)cx / (GLfloat)cy;

//	psizex = (GLfloat)cx;
//	psizey = (GLfloat)cy;

	///here set the width// KBH!!!!!!!


	glMatrixMode(GL_PROJECTION); //원근 투영 임을 선언

	glLoadIdentity();

	gluPerspective(60.0f, fAspect, 1.0f, 10000.0f); //관측 공간을 설정 해 준다
	//gluPerspective(60.0f, fAspect, 1.0f, 40.0f); //관측 공간을 설정 해 준다

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}

int CMFCGLView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	
	m_pDC = new CClientDC(this); // Get device context
	
	if(NULL==m_pDC)
	{
//		AfxMessageBox("Cannot get DC.\n");
		return FALSE;
	}
	
	if(!SetupPixelFormat())
	{  // Set up the fixel format
//		::AfxMessageBox("SetupPixelFormat fail.\n");
		return FALSE;
	}

	// Get rendering context
	if(0==(m_hRC=wglCreateContext(m_pDC->GetSafeHdc())))
	{
//		::AfxMessageBox("wglCreateContext fail.\n");
		return FALSE;
	}
	
	// make rendering context
	if(FALSE==wglMakeCurrent(m_pDC->GetSafeHdc(), m_hRC))
	{
//		::AfxMessageBox("wglMakeCurrent fail.\n");
		return FALSE;
	}
	return 0;
}

BOOL CMFCGLView::SetupPixelFormat(PIXELFORMATDESCRIPTOR * pPFD)
{
	PIXELFORMATDESCRIPTOR pfd =
	{
		sizeof(PIXELFORMATDESCRIPTOR),  // PIXELFORMATDESCRIPTOR SIZE
		1,                              // Version
		PFD_DRAW_TO_WINDOW |            // support window
		PFD_SUPPORT_OPENGL
		| PFD_DOUBLEBUFFER,             // support OpenGL
		PFD_TYPE_RGBA,                  // RGBA type
		24,                             // 24-bit color depth
		0, 0, 0, 0, 0, 0,               // color bits ignored
		0,                              // no alpha buffer
		0,                              // shift bit ingored
		0,                              // no accumulation buffer
		0, 0, 0, 0,                     // accum bits ingored
		16,                             // 16-bit z-buffer
		0,                              // no stencil buffer
		0,                              // no auxiliary buffer
		PFD_MAIN_PLANE,                 // main layer
		0,                              // reserved
		0, 0, 0                         // layer masks ignored
	};

    int pixelformat;

	PIXELFORMATDESCRIPTOR* pPFDtoUse;


	// let the user override the default pixel format
	pPFDtoUse = (0 == pPFD)? &pfd : pPFD;

	if(0==(pixelformat=::ChoosePixelFormat(m_pDC->GetSafeHdc(), pPFDtoUse))) 
	{
//		::AfxMessageBox("ChoosePixelFormat failed.");
		return FALSE;
	}
	
	if(FALSE==::SetPixelFormat(m_pDC->GetSafeHdc(), pixelformat, pPFDtoUse)) 
	{
//		::AfxMessageBox("SetPixelFormat failed.");
		return FALSE;
	}

	return TRUE;
}

GLvoid CMFCGLView::BuildFont(GLvoid)
{
    HFONT font;                           // <1>
    base = glGenLists(24);                // <2>
    font = CreateFont(-12,                // <3-1>  // fontsize
                        0,
                        0,
                        0,
                  FW_BOLD,                // <3-2>
                    FALSE,                // <3-3>
                    FALSE,                // <3-4>
                    FALSE,                // <3-5>
             ANSI_CHARSET,                // <3-6>
            OUT_TT_PRECIS,
      CLIP_DEFAULT_PRECIS,
      ANTIALIASED_QUALITY,
FF_DONTCARE|DEFAULT_PITCH,
                          "Courier New"); // <3-6>
   
    SelectObject(m_pDC->GetSafeHdc(), font);              // <4>
    wglUseFontBitmaps(m_pDC->GetSafeHdc(), 32, 96, base); // <5>
}
GLvoid CMFCGLView::KillFont(GLvoid)
{
    glDeleteLists(base, 96);
}

GLvoid CMFCGLView::glPrint(const char *text)
{
    glPushAttrib(GL_LIST_BIT);                         //<1>
    glListBase(base - 32);                             //<2>
    glCallLists(strlen(text), GL_UNSIGNED_BYTE, text); //<3>
    glPopAttrib();                                     //<4>
}

int CMFCGLView::InitGL(GLvoid)
{
    glShadeModel(GL_SMOOTH);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
   
    ///////////////////////// NEW //////////////////////////
  //  BuildFont();  //this exists in the showwindows function
    ///////////////////////// NEW //////////////////////////
   
    return TRUE;
}


void CMFCGLView::DrawScene(void)
{


	
	// LIGHT 
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// DEPTH
	glEnable(GL_DEPTH_TEST);
	
//	glClearColor(1.0, 1.0, 1.0, 0.0);

	glClearColor(0.8, 0.8, 0.8, 0.0);


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// MATERIAL PROPERTY
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);

	// MATRIX MODE
	glMatrixMode(GL_MODELVIEW);

	//COORDINATE ROTATION
	glLoadIdentity();
 
	gluLookAt(0.0f, 0.0f, -4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);   //카메라 설정
	//카메라의 위치는 (0, 0, -1000) 에 있으며, 카메라가 바라보는 곳은 (0, 10, 0) 이다.
	//카메라의 Up Vector 는 (0, 1, 0) 이다.


	glColor3f(0.0f, 0.0f, 1.0f);                            //버텍스의 색깔을 흰색으로 셋팅


	glScalef(-1.0f, 1.0f, 1.0f); //Return origin

	glTranslatef(float(m_xTranslate),float(m_yTranslate),float(m_zTranslate));
	glRotatef(float(m_xRotate), 1.0, 0.0, 0.0);
	glRotatef(float(m_yRotate), 0.0, 1.0, 0.0);
    glRotatef(float(m_zRotate), 0.0, 0.0, 1.0);

/*
	 glTranslatef(0.0f, -0.3f, 0.0f); //  glColor3f(0.0f, 1.0f, 0.0f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("OpenGL Template by Bohung Kim, University of Ulsan");  //<3>
 	 glTranslatef(-0.0f, 0.3f, 0.0f);


	glTranslatef(0.0f, 0.0f, 0.0); //Target point
	auxSolidSphere(1.0f);         //Draw sphere
	glTranslatef(0.0f, 0.0f, 0.0);  //Return origin

	glTranslatef(2.5f, 0.0f, 0.0); //Target point
	auxSolidSphere(1.0f);         //Draw sphere
	glTranslatef(-2.5f, 0.0f, 0.0);  //Return origin

	glTranslatef(0.0f, 2.5f, 0.0); //Target point
	auxSolidSphere(1.0f);         //Draw sphere
	glTranslatef(0.0f, -2.5f, 0.0); //Return origin

	glTranslatef(0.0f, 0.0f, 2.5f); //Target point
	auxSolidSphere(1.0f);         //Draw sphere
	glTranslatef(0.0f, 0.0f, -2.5f); //Return origin


	
		
	glTranslatef(glcrossx, glcrossy,0.0f); //Target point
	auxSolidSphere(1.0f);         //Draw sphere
	glTranslatef(-glcrossx, -glcrossy,0.0f); //Return origin







	glColor3f(1.0f, 1.0f, 1.0f); //Set the sphere color to RED

	glBegin(GL_LINES);

		 glVertex3f(0.0f, 0.0f, 0.0f);
		 glVertex3f(float(m_fW), 0.0f,0.0f);
 		 glVertex3f(float(m_fW), 0.0f,0.0f);
		 glVertex3f(float(m_fW), float(m_fH),0.0f);
		 glVertex3f(float(m_fW), float(m_fH),0.0f);
		 glVertex3f(0.0f, float(m_fH),0.0f);
		 glVertex3f(0.0f, float(m_fH),0.0f);
		 glVertex3f(0.0f, 0.0f, 0.0f);

	glColor3f(1.25f, 1.25f, 1.25f);
		 glVertex3f(0.0f, 0.0f, float(m_fZ));
		 glVertex3f(float(m_fW), 0.0f, float(m_fZ));
 		 glVertex3f(float(m_fW), 0.0f,float(m_fZ));
		 glVertex3f(float(m_fW), float(m_fH),float(m_fZ));
		 glVertex3f(float(m_fW), float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, 0.0f, float(m_fZ));

	glColor3f(1.65f, 1.65f, 1.65f);
		 glVertex3f(0.0f, 0.0f, 0.0f);
		 glVertex3f(0.0f, 0.0f, float(m_fZ));
 		 glVertex3f(float(m_fW), 0.0f,0.0f);
		 glVertex3f(float(m_fW), 0.0f,float(m_fZ));
		 glVertex3f(float(m_fW), float(m_fH),0.0f);
		 glVertex3f(float(m_fW), float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, float(m_fH),0.0f);
		 glVertex3f(0.0f, float(m_fH), float(m_fZ));

	 glEnd();


*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
/////////////////////////////////////////////////////////////// Based on CMolecular Object

	
	// distance  : THIS IS FOR THE UNIT DISTANCE OF 1.0 IN THIS PROGRAM
	glTranslatef(-1.5f, 0.25f, 0.0); //Target point
	auxSolidSphere(0.1f);         //Draw sphere
	glTranslatef(1.5f, -0.25f, 0.0);  //Return origin

	glTranslatef(-1.5f, 1.25f, 0.0); //Target point
	auxSolidSphere(0.1f);         //Draw sphere
	glTranslatef(1.5f, -1.25f, 0.0); //Return origin

	glColor3f(0.75f, 0.0, 0.0);
	//Draw velocity marker //Left BC
	//for (int i=0; i<(int(m_fH)+1) ; i++)
	/*
	for (int i=0; i<(int(m_fH/m_ygsize)) ; i++)
	{
	glTranslatef(float(m_aMVel[i].m_fX),float(m_aMVel[i].m_fY), 0.0); //Target point
	auxSolidSphere(0.2f);         //Draw sphere
	glTranslatef(float(-m_aMVel[i].m_fX),float(-m_aMVel[i].m_fY), 0.0); //Return origin
	}*/
	
/*
	//Draw Right BC
	for (int i=0; i<int(m_fH/m_ygsize) ; i++)
	{
	glTranslatef(10.0f,float(m_aMVel[i].m_fY), 0.0); //Target point
	auxSolidSphere(0.1f);         //Draw sphere
	glTranslatef(-10.0f,float(-m_aMVel[i].m_fY), 0.0); //Return origin
	}
*/
///////////////////////////////////////////////////////////////////////////////////

	//Draw all molecules

	glColor3f(0.0, 0.0, float(m_dautocolor)); //Set the sphere color to BLUE
	//Drawing top wall
	/*
	for(int ii = 0; ii<m_nwMs ; ii++)
	{
		//if (abs(int(m_aTwM[ii].m_fay*1000000000000000)) != 0.0f)
		//{glColor3f(0.75f, 0.75f, 0.0);}
		//else
		//{glColor3f(0.0, 0.0, float(m_dautocolor));}
		
		 if (abs(int(m_aTwM[ii].m_fay*10)) != 0.0f)
		{glColor3f(0.75f, 0.75f, 0.0);}
		else
		{glColor3f(0.0, 0.0, float(m_dautocolor));}


		glTranslatef(float(m_aTwM[ii].m_fX),float(m_aTwM[ii].m_fY),float(m_aTwM[ii].m_fZ));   //Target point
		if (m_nWire == 0) {auxSolidSphere(m_fMrad);}
		else{auxWireSphere(m_fMrad);}
		glTranslatef(float(-m_aTwM[ii].m_fX),float(-m_aTwM[ii].m_fY),float(-m_aTwM[ii].m_fZ)); //Return to origin

	//Drawing bottom wall
	 if (abs(int(m_aBwM[ii].m_fay*10)) != 0.0f)
		{glColor3f(0.75f, 0.75f, 0.0);}
		else
		{glColor3f(0.0, 0.0, float(m_dautocolor));}


		glTranslatef(float(m_aBwM[ii].m_fX),float(m_aBwM[ii].m_fY),float(m_aBwM[ii].m_fZ));    //Target point
		if (m_nWire == 0) {auxSolidSphere(m_fMrad);}
		else{auxWireSphere(m_fMrad);}
		glTranslatef(float(-m_aBwM[ii].m_fX),float(-m_aBwM[ii].m_fY),float(-m_aBwM[ii].m_fZ)); //Return to origin
	}

	glColor3f(0.0, 0.75f, 0.0); // Change sphere color to GREEN

	//Drawing inner Molecule
	for( ii = 0; ii<m_nMs ; ii++)
	{
		if (ii > (m_nMs/2-1)) glColor3f(0.75f, 0.75f, 0.75f);
		glTranslatef(float(m_aM[ii].m_fX),float(m_aM[ii].m_fY),float(m_aM[ii].m_fZ)); //target point
		if (m_nWire == 0) {auxSolidSphere(m_fMrad);}
		else{auxWireSphere(m_fMrad);}
		glTranslatef(float(-m_aM[ii].m_fX),float(-m_aM[ii].m_fY),float(-m_aM[ii].m_fZ)); //return to origin
	}
*/

	int xx,yy,zz,nn;
	double dNCX,dNCY,dNCZ;
	dNCX = m_fW/m_cellsize;
	dNCY = m_fH/m_cellsize;
	dNCZ = m_fZ/m_cellsize;

	for (xx = 0; xx < dNCX + 2 ; xx++)
	{
		for ( yy = 0; yy < dNCY + 2 ; yy++)
		{
			for ( zz = 0; zz < dNCZ + 2 ; zz++)
			{
				for (nn = 0; nn< Cell[xx][yy][zz].NM ; nn++) 
				{
				glTranslatef(float(Cell[xx][yy][zz].x[nn]),float(Cell[xx][yy][zz].y[nn]),float(Cell[xx][yy][zz].z[nn])); //target point
				if (Cell[xx][yy][zz].mtype[nn]== 1) //fluid
				{
					glColor3f(0.0, 0.75f, 0.0); // Change sphere color to GREEN
					auxSolidSphere(m_fMrad);
				}

				if (Cell[xx][yy][zz].mtype[nn] == 4) // thrmal wall
				{
					if (m_nSp != 0)
					{glColor3f(1.00f, 0.65f, 0.0f);}
					else
					{glColor3f(0.75f, 0.75f, 0.0f);}
					auxSolidSphere(m_fMrad);
				}


				if (Cell[xx][yy][zz].mtype[nn] == 2) // wall
				{
					if (m_nWire != 1) 
					{
						if (m_nWire != 3) 
						{


						if (m_nSp != 0)
						{glColor3f(0.15f, 0.25f, 0.9f);}
						else
						{glColor3f(0.0, 0.0f, 0.75f);}

					auxSolidSphere(m_fMrad);
						}
					}
				}
				if (Cell[xx][yy][zz].mtype[nn] == 3) //image
				{
					if (m_nWire != 2) 
					{
						if (m_nWire != 3) 
						{
					glColor3f(0.0, 0.2f, 0.0f); // Change sphere color to GREEN
					//auxWireSphere(m_fMrad/2);
					auxWireSphere(m_fMrad);
					//	auxSolidSphere(m_fMrad-0.02);
						}
					}
				}
				glTranslatef(float(-Cell[xx][yy][zz].x[nn]),float(-Cell[xx][yy][zz].y[nn]),float(-Cell[xx][yy][zz].z[nn])); //return to origin
				}

			}
		}
	}


//str.Format(_T("Floating point: %.2f\n"), 12345.12345);
//_tprintf(_T("%s"), (LPCTSTR) str);

//str.Format(_T("Left-justified integer: %.6d\n"), 35);
//_tprintf(_T("%s"), (LPCTSTR) str);

//str.Format(IDS_SCORE, 5, 3);
//_tprintf(_T("%s"), (LPCTSTR) str);


//Floating point: 12345.12
//Left-justified integer: 000035
//Penguins: 5
//Flyers  : 3



	//m_T = temperature

	//char buffer1[20];
	//_gcvt( m_T, 3, buffer1 );
	//m_strTemp = buffer1;
	m_strTemp.Format("%.2f",m_T);
	m_strTemp = m_strTemp + "K";


	//char buffer[20];
	//_itoa(m_nMs, buffer, 10 );
	//char buffer2[20];
	//_gcvt( m_nMs, 3, buffer2 );
	//_itoa(m_nMs, buffer2, 10 );
	//m_strNM = buffer2;
	m_strNM.Format("%d",m_nMs);
	m_strNM = "M #: " + m_strNM;

	//char buffer3[20];
	//_gcvt( d_Twalltemp, 3, buffer3 );
	//m_strTWallTemp = buffer3;
	m_strTWallTemp.Format("%.2f",d_Twalltemp);
	m_strTWallTemp = "Top Wall Temp: " + m_strTWallTemp;

	//char buffer4[20];
	//_gcvt( d_Bwalltemp, 3, buffer4 );
	//m_strBWallTemp = buffer4;
	m_strBWallTemp.Format("%.2f",d_Bwalltemp);
	m_strBWallTemp = "Bottom Wall Temp: " + m_strBWallTemp;

	//char bufferP[20];
	//_gcvt( m_P, 3, bufferP );
	//m_strP = bufferP;
	m_strP.Format("%.2f",m_P);
	m_strP = "Pressure(KPA): " + m_strP;


	//CString NumberofStep,totalLoop,heatfluxx,heatfluxy,heatfluxz;

	//char NumberofStep[20];
	//_gcvt( d_ntimestep, 6, NumberofStep );
m_str_NumberofStep.Format("%.0f",d_ntimestep);
	 glTranslatef(16.5f, 10.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_str_NumberofStep);  //<3>
 	 glTranslatef(-16.5f, -10.1f, 0.0f);


  	//char totalLoop[20];
 	//_gcvt( d_totalLoop, 6, totalLoop );
m_str_totalLoop.Format("%.0f",d_totalLoop);
	 glTranslatef(float(m_fW) + 2.5f, 11.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_str_totalLoop);  //<3>
 	 glTranslatef(-float(m_fW) + -2.5f, -11.1f, 0.0f);

	 //char heatfluxx[20];
	//_gcvt( m_heatfluxx, 6, heatfluxx );
m_str_heatfluxx.Format("%.2f",m_heatfluxx);
	 glTranslatef(float(m_fW) + 3.5f, 11.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_str_heatfluxx);  //<3>
 	 glTranslatef(-float(m_fW) - 3.5f, -11.1f, 0.0f);

	 //char heatfluxy[20];
	//_gcvt( m_heatfluxy, 6, heatfluxy );
m_str_heatfluxy.Format("%.2f",m_heatfluxy);
	 glTranslatef(float(m_fW) + 3.5f, 10.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_str_heatfluxy);  //<3>
 	 glTranslatef(-float(m_fW) - 3.5f, -10.6f, 0.0f);

	 //char heatfluxz[20];
	//_gcvt( m_heatfluxz, 6, heatfluxz );
m_str_heatfluxz.Format("%.2f",m_heatfluxz);
	 glTranslatef(float(m_fW) + 3.5f, 10.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_str_heatfluxz);  //<3>
 	 glTranslatef(-float(m_fW) - 3.5f, -10.1f, 0.0f);




	if (m_time == 0) 
	{
     glTranslatef(12.5f, 5.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strBWallTemp);  //<3>
 	 glTranslatef(-12.5f, -5.6f, 0.0f);

      glTranslatef(14.5f, 7.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strP);  //<3>
 	 glTranslatef(-14.5f, -7.6f, 0.0f);



     glTranslatef(12.5f, 6.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strTWallTemp);  //<3>
 	 glTranslatef(-12.5f, -6.1f, 0.0f);

     glTranslatef(12.5f, 6.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("Simulation Off");  //<3>
 	 glTranslatef(-12.5f, -6.6f, 0.0f);
	 
	 glTranslatef(12.5f, 7.1f, -0.0f);   glColor3f(1.75f, 1.0f, 0.0f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strTemp);  //<3>
 	 glTranslatef(-12.5f, -7.1f, 0.0f);

	 glTranslatef(12.5f, 7.6f, -0.0f);   glColor3f(1.75f, 1.0f, 0.0f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strNM);  //<3>
 	 glTranslatef(-12.5f, -7.6f, 0.0f);
	}else
	{
     glTranslatef(12.5f, 5.6f, -0.0f);   glColor3f(1.75f, 2.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strBWallTemp);  //<3>
 	 glTranslatef(-12.5f, -5.6f, 0.0f);

	 glTranslatef(14.5f, 7.6f, -0.0f);   glColor3f(1.75f, 2.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strP);  //<3>
 	 glTranslatef(-14.5f, -7.6f, 0.0f);


     glTranslatef(12.5f, 6.1f, -0.0f);   glColor3f(1.75f, 2.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strTWallTemp);  //<3>
 	 glTranslatef(-12.5f, -6.1f, 0.0f);

     glTranslatef(12.5f, 6.6f, -0.0f);   glColor3f(1.75f, 2.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("Simulation On");  //<3>
 	 glTranslatef(-12.5f, -6.6f, 0.0f);

	 glTranslatef(12.5f, 7.1f, -0.0f);   glColor3f(2.75f, 0.2f, 0.0f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strTemp);  //<3>
 	 glTranslatef(-12.5f, -7.1f, 0.0f);

	 glTranslatef(12.5f, 7.6f, -0.0f);   glColor3f(2.75f, 0.2f, 0.0f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(m_strNM);  //<3>
 	 glTranslatef(-12.5f, -7.6f, 0.0f);
	}
	//Temperature

	//CString BT1,BT2,BT3,TT1,TT2,TT3;
	//char buffer5[20];char buffer6[20];char buffer7[20];
	//char buffer8[20];char buffer9[20];char buffer10[20];
	//_gcvt( d_Bwalltemp1, 3, buffer5 );
	//_gcvt( d_Bwalltemp2, 3, buffer6 );
	//_gcvt( d_Bwalltemp3, 3, buffer7 );
	//_gcvt( d_Twalltemp1, 3, buffer8 );
	//_gcvt( d_Twalltemp2, 3, buffer9 );
	//_gcvt( d_Twalltemp3, 3, buffer10 );

	BT1.Format("%.2f",d_Bwalltemp1);
	BT2.Format("%.2f",d_Bwalltemp2);
	BT3.Format("%.2f",d_Bwalltemp3);
	TT1.Format("%.2f",d_Twalltemp1);
	TT2.Format("%.2f",d_Twalltemp2);
	TT3.Format("%.2f",d_Twalltemp3);

    glTranslatef(15.5f, 3.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
    glRasterPos2f(0.0f, 0.0f);   glPrint(TT3);  //<3>
 	glTranslatef(-15.5f, -3.6f, 0.0f);
	glTranslatef(15.5f, 3.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
    glRasterPos2f(0.0f, 0.0f);   glPrint(TT2);  //<3>
 	glTranslatef(-15.5f, -3.1f, 0.0f);
	glTranslatef(15.5f, 2.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
    glRasterPos2f(0.0f, 0.0f);   glPrint(TT1);  //<3>
 	glTranslatef(-15.5f, -2.6f, 0.0f);
	glTranslatef(15.5f, 2.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
    glRasterPos2f(0.0f, 0.0f);   glPrint(BT1);  //<3>
 	glTranslatef(-15.5f, -2.1f, 0.0f);
	glTranslatef(15.5f, 1.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
    glRasterPos2f(0.0f, 0.0f);   glPrint(BT2);  //<3>
 	glTranslatef(-15.5f, -1.6f, 0.0f);
	glTranslatef(15.5f, 1.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
    glRasterPos2f(0.0f, 0.0f);   glPrint(BT3);  //<3>
 	glTranslatef(-15.5f, -1.1f, 0.0f);





//	if (m_time == 0) 
//	{

	int ii;

//	char tbuffer[10][20];
//	char dbuffer[10][20];

	for (ii = 0; ii< 10; ii++)
	{
		if (dG[ii] < 10000.0)
		{
		//_gcvt( dG[ii], 3, tbuffer[ii] );
		//strG[ii] = tbuffer[ii];
		strG[ii].Format("%.2f",dG[ii]);

		//_gcvt( denCount[ii], 3, dbuffer[ii] );
		//strGDen[ii] = dbuffer[ii];
		strGDen[ii].Format("%.0f",denCount[ii]);
		}
	}

	for (ii = 0; ii< 10; ii++)
	{
    glTranslatef(float(m_fW) + 2.5f, 0.1f+float(ii)*0.5f, 0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);  
	 glPrint(strG[ii]);  //<3>
	 //glPrint("d");  //<3>
 	 glTranslatef(-float(m_fW) + -2.5f, -(0.1f+float(ii)*0.5f), 0.0f);
	}

	for (ii = 0; ii< 10; ii++)
	{
    glTranslatef(float(m_fW) + 4.0f, 0.1f+float(ii)*0.5f, 0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);  
	 glPrint(strGDen[ii]);  //<3>
	 //glPrint("d");  //<3>
 	 glTranslatef(-float(m_fW) + -4.0f, -(0.1f+float(ii)*0.5f), 0.0f);
	}


     glTranslatef(12.5f, 10.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("iCPU, sec/loop");  //<3>
 	 glTranslatef(-12.5f, -10.6f, 0.0f);


//	 char stricpu[20];
//	 _itoa(iCPU, stricpu, 10 );

// 	 char strsecperloop[20];
//	 _gcvt( duration , 6, strsecperloop );

	 CString strsecperloop,stricpu;

	strsecperloop.Format("%.3f",duration);
	stricpu.Format("%d",iCPU);

	 glTranslatef(15.5f, 10.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(stricpu);  //<3>
 	 glTranslatef(-15.5f, -10.6f, 0.0f);

	 glTranslatef(16.5f, 10.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint(strsecperloop);  //<3>
 	 glTranslatef(-16.5f, -10.6f, 0.0f);


     glTranslatef(12.5f, 10.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("Argon [Ar]");  //<3>
 	 glTranslatef(-12.5f, -10.1f, 0.0f);

	 glTranslatef(12.5f, 9.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("A# 18");  //<3>
 	 glTranslatef(-12.5f, -9.6f, 0.0f);

	 glTranslatef(12.5f, 9.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("AM 39.948");  //<3>
 	 glTranslatef(-12.5f, -9.1f, 0.0f);

	 glTranslatef(12.5f, 8.6f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("MT 84K ");  //<3>
 	 glTranslatef(-12.5f, -8.6f, 0.0f);

	 glTranslatef(12.5f, 8.1f, -0.0f);   glColor3f(1.75f, 1.75f, 1.75f);      //<1>
     glRasterPos2f(0.0f, 0.0f);   glPrint("BT 89K ");  //<3>
 	 glTranslatef(-12.5f, -8.1f, 0.0f);




//	}

	 //Draw TEmp MArker

		 	glColor3f(0.75f, 0.0f, 0.0);
				glTranslatef(11.0f+float(m_Gcount)/50.0f, float(m_Tsc)+4.5f, 0.0f); //Target point
				auxSolidSphere(0.08f);         //Draw sphere
				glTranslatef(-11.0f-float(m_Gcount)/50.0f, -float(m_Tsc)-4.5f, 0.0f);  //Return origin
				glColor3f(5.75f, 0.0f, 0.0);


	 glBegin(GL_LINE_STRIP);
		 for (ii=1; ii<m_Gcount ;ii++)
		 {
	

			 glVertex3f(11.0f+float(ii)/50.0f,float(m_Tsca[ii])+4.5f, 0.0f); 
		 }
	 glEnd();


	//Draw Kin MArker
	glColor3f(0.0f, 0.0f, 0.75f);
	glTranslatef(11.0f+float(m_Gcount)/50.0f, float(m_Kinsc)+1.0f, 0.0f); //Target point
	auxSolidSphere(0.08f);         //Draw sphere
	glTranslatef(-11.0f-float(m_Gcount)/50.0f, -float(m_Kinsc)-1.0f, 0.0f);  //Return origin
	glColor3f(0.0f, 0.0f, 5.75f);
	 glBegin(GL_LINE_STRIP);
		 for (ii=1; ii<m_Gcount ;ii++)
		 { glVertex3f(11.0f+float(ii)/50.0f,float(m_Kinsca[ii])+1.01f, 0.0f); }
	 glEnd();


	//Draw Pot MArker
	glColor3f(0.0f, 0.0f, 0.0f);
	glTranslatef(11.0f+float(m_Gcount)/50.0f, float(m_rf*m_Potsc)+1.0f, 0.0f); //Target point
	auxSolidSphere(0.08f);         //Draw sphere
	glTranslatef(-11.0f-float(m_Gcount)/50.0f,-float(m_rf*m_Potsc)-1.0f, 0.0f);  //Return origin
	glColor3f(0.0f, 0.0f, 0.0f);
	 glBegin(GL_LINE_STRIP);
		 for (ii=1; ii<m_Gcount ;ii++)
		 { glVertex3f(11.0f+float(ii)/50.0f,float(m_rf*m_Potsca[ii])+1.0f, 0.0f); }
	 glEnd();



	//Draw Tot MArker
	glColor3f(0.0f, 1.0f, 0.0f);
	glTranslatef(11.0f+float(m_Gcount)/50.0f, float(m_rf*m_Totsc)+1.0f, 0.0f); //Target point
	auxSolidSphere(0.08f);         //Draw sphere
	glTranslatef(-11.0f-float(m_Gcount)/50.0f, -float(m_rf*m_Totsc)-1.0f, 0.0f);  //Return origin
	glColor3f(0.0f, 1.0, 0.0f);
	 glBegin(GL_LINE_STRIP);

		 for (ii=1; ii<m_Gcount ;ii++)
		 { glVertex3f(11.0f+float(ii)/50.0f,float(m_rf*m_Totsca[ii])+1.0f, 0.0f); }
	 glEnd();


	 // Draw Plot
	 glColor3f(2.75f, 2.75f, 2.75f);
     //X
	 glBegin(GL_LINES);
	 glVertex3f(11.0f, 1.0f, 0.0f);
	 glVertex3f(19.5f, 1.0f, 0.0f);
	 glEnd();
	//Temp
 	 glBegin(GL_LINES);
	 glVertex3f(11.0f, 4.5f, 0.0f);
	 glVertex3f(19.5f, 4.5f, 0.0f);
	 glEnd();

	 //Y
	 glBegin(GL_LINES);
	 glVertex3f(11.0f, -2.0f, 0.0f);
	 glVertex3f(11.0f, 12.0f, 0.0f);
	 glEnd();

	 //grid

	 glColor3f(0.75f, 0.75f, 0.0f);
	 for (int kk = 0; kk < 27 ; kk++)
	 {
	 float grid;
	 grid = float(kk)*0.5f-1.5f;
	 glBegin(GL_LINES);
	 glVertex3f(11.0f, grid, 0.0f);
	 glVertex3f(19.5f, grid, 0.0f);
	 glEnd();
	 }

  /*
	 
*/
	glColor3f(2.75f, 2.75f, 2.75f); //Set the sphere color to RED

	glBegin(GL_LINES);

		 glVertex3f(0.0f, 0.0f, 0.0f);
		 glVertex3f(float(m_fW), 0.0f,0.0f);
 		 glVertex3f(float(m_fW), 0.0f,0.0f);
		 glVertex3f(float(m_fW), float(m_fH),0.0f);
		 glVertex3f(float(m_fW), float(m_fH),0.0f);
		 glVertex3f(0.0f, float(m_fH),0.0f);
		 glVertex3f(0.0f, float(m_fH),0.0f);
		 glVertex3f(0.0f, 0.0f, 0.0f);

	glColor3f(1.25f, 1.25f, 1.25f);
		 glVertex3f(0.0f, 0.0f, float(m_fZ));
		 glVertex3f(float(m_fW), 0.0f, float(m_fZ));
 		 glVertex3f(float(m_fW), 0.0f,float(m_fZ));
		 glVertex3f(float(m_fW), float(m_fH),float(m_fZ));
		 glVertex3f(float(m_fW), float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, 0.0f, float(m_fZ));

	glColor3f(1.65f, 1.65f, 1.65f);
		 glVertex3f(0.0f, 0.0f, 0.0f);
		 glVertex3f(0.0f, 0.0f, float(m_fZ));
 		 glVertex3f(float(m_fW), 0.0f,0.0f);
		 glVertex3f(float(m_fW), 0.0f,float(m_fZ));
		 glVertex3f(float(m_fW), float(m_fH),0.0f);
		 glVertex3f(float(m_fW), float(m_fH),float(m_fZ));
		 glVertex3f(0.0f, float(m_fH),0.0f);
		 glVertex3f(0.0f, float(m_fH), float(m_fZ));

	 glEnd();

/*
	 for (int bb = 0; bb < m_fW/m_cellsize+1 ; bb++)
	 {
		 for (int aa = 0; aa < m_fW/m_cellsize+1 ; aa++)
		 {
		 glBegin(GL_LINES);
		 glVertex3f(float(m_cellsize*float(aa)), 0.0f, -float(m_cellsize*float(bb)));
		 glVertex3f(float(m_cellsize*float(aa)), float(m_fH), -float(m_cellsize*float(bb)));
		 glEnd();
		 }

		 for ( aa = 0; aa < m_fH/m_cellsize+1 ; aa++)
		 {
		 glBegin(GL_LINES);
		 glVertex3f(0.0f, float(m_cellsize*float(aa)), -float(m_cellsize*float(bb)));
		 glVertex3f(float(m_fW),  float(m_cellsize*float(aa)), -float(m_cellsize*float(bb)));
		 glEnd();
		 }
	 }
*/

/*
주기율표 제0족에 속하는 비활성기체원소. 
 
   원소기호  Ar 
   원자번호  18 
   원자량  39.948 
   녹는점  -189.2℃ 84K
   끓는점  -185.7℃ 89K
   비중  1.7834g/ℓ 
 
 본문 

 1894년에 영국의 J.W.레일리가 W.램지의 협력을 얻어 공기에서
 산소를 제거하고 얻은 질소의 비중과, 질소화합물을 분해하여
 얻은 질소의 비중이 다른 점에 착안하여, 공기에서 이 물질을
 분리시켜 발견하였다.  이 물질은 그때까지 알려진 원소와 달리
 화학적 성질이 극히 활발하지 않아 모든 물질과 반응하지 않았던
 데서, 그리스어인 'an ergon(게으름쟁이)'을 따서 명명되었다.
 */




	/////////////////////////////////////////

	





	glFlush();



}

void CMFCGLView::OnDestroy()
{

	CView::OnDestroy();
	
	// TODO: Add your message handler code here

	if(FALSE==wglMakeCurrent(0, 0)) {  //Release current rendering context
		//::AfxMessageBox("wglMakeCurrent Fail.\n");
	}
	
	if(m_hRC && (FALSE==wglDeleteContext(m_hRC))) { // Removing rendering context
		//::AfxMessageBox("wglDeleteContext Fail.\n");
	}

	if(m_pDC) delete m_pDC;  // Remove the device context


	// TODO: Add your message handler code here
}

void CMFCGLView::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

//	nRotate =1;  //0=fix, 1=Rotate, 
//	m_xRotate = m_xRotate+300.0;
//	m_xRotate=(double)point.y*30000.0;
//	m_yRotate=(double)point.x*30000.0;
//
//	DrawScene();
//	SwapBuffers(m_hDC);

	CView::OnLButtonDown(nFlags, point);
}

void CMFCGLView::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

//	nRotate =0;  //0=fix, 1=Rotate, 
	CView::OnLButtonUp(nFlags, point);
}

void CMFCGLView::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

//	if (nRotate==1)
//	{

//	m_xRotate = m_xRotate+30.0;
//	m_xRotate=(double)point.y;
//	m_yRotate=(double)point.x;
//

	
	int    vp[4];
	double mv[16], pm[16];
	float  winX, winY, winZ;
	double glX, glY, glZ;
	CPoint cursor;

	GetCursorPos(&cursor);
	ScreenToClient(&cursor);

	glGetIntegerv(GL_VIEWPORT, vp);
	glGetDoublev(GL_MODELVIEW_MATRIX, mv);
	glGetDoublev(GL_PROJECTION_MATRIX, pm);

	winX = cursor.x;
	winY = cursor.y;

	glReadPixels(winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);

	gluUnProject(winX, winY, 0.0f, mv, pm, vp, &glX, &glY, &glZ);
	 
	glcrossx = 15.0f*(glX-8.000);
	glcrossy = 15.0f*(-glY+4.7f);
	//for future use	

		
	if ( nFlags == MK_LBUTTON )
	{
	m_xRotate=(double)point.y;
    m_yRotate=(double)point.x;
	}

	if ( nFlags == MK_RBUTTON)
	{
	m_xTranslate= m_xTranslate + ((double)point.x-500)/3500.0 ;
	m_yTranslate= m_yTranslate - ((double)point.y-300)/3500.0 ;
	}





	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
	//SwapBuffers(m_hDC);

	//}
	

//	CView::OnMouseMove(nFlags, point);
}

void CMFCGLView::OnRButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

	CView::OnRButtonDown(nFlags, point);
}

void CMFCGLView::OnTimer(UINT_PTR nIDEvent)
{
	// TODO: Add your message handler code here and/or call default
	//OMPTest();

		// TODO: Add your message handler code here and/or call default
	if (m_time == 1) 
	{
		OnRun();
	//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}
	CView::OnTimer(nIDEvent);

}

void CMFCGLView::OMPTest()
{

//	clock_t start, finish;
//	double  duration;

	// Get the number of processors in this system
	//iCPU = omp_get_num_procs();

	//iCPU = 1;
	// Now set the number of threads
	omp_set_num_threads(iCPU);
	
//	start = clock();

		#pragma omp parallel for
		for(int y = 0; y < 200000; y++) 
			{ 
				for(int x = 0; x< 200000; x++) 
				{ 
						 x = x;
				} 
			}

	m_xRotate=m_xRotate+10.0;
	m_yRotate=m_yRotate+10.0;

	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());		//SwapBuffers(m_hDC);


//	finish = clock();

//	duration = (double)(finish - start) / CLOCKS_PER_SEC;

}

void CMFCGLView::OnAutoRun()
{

	OnRun();
	// TODO: Add your command handler code here
	//SetTimer(1,1,NULL); //SetTimer(timer ID,milisecond,NULL)
}

void CMFCGLView::OnAutoStop()
{
	// TODO: Add your command handler code here
	//KillTimer(1);
}

void CMFCGLView::OnAutoSinglecpu()
{
	// TODO: Add your command handler code here
	iCPU = 1;

	omp_set_num_threads(iCPU);

}

void CMFCGLView::OnAutoAddcpu()
{
	// TODO: Add your command handler code here
	int max = omp_get_num_procs();
	iCPU = iCPU + 1;
	if (iCPU == max) iCPU =max;

	omp_set_num_threads(iCPU);
}

void CMFCGLView::OnAutoMaxcpu()
{
	// TODO: Add your command handler code here
	iCPU = omp_get_num_procs();

	omp_set_num_threads(iCPU);

}

void CMFCGLView::OnAutoRemovecpu()
{
	// TODO: Add your command handler code here
	iCPU = iCPU - 1;
	if (iCPU ==0 ) iCPU =1;

	omp_set_num_threads(iCPU);

}

void CMFCGLView::OnShowWindow(BOOL bShow, UINT nStatus)
{
	BuildFont();
	CView::OnShowWindow(bShow, nStatus);
}


void CMFCGLView::TmpFileSave()
{

	m_file.Open("MDtmpSave.txt",CFile::modeCreate | CFile::modeWrite);

	CString temp;
//	char buffer[20];

	//_itoa(m_nMs, buffer, 10 );
	//temp = CString(buffer);
	temp.Format("%d",m_nMs);

	m_file.WriteString(temp+"\n");  //Write Number of Moles

	//char buffer1[20],buffer2[20],buffer3[20],buffer4[20],buffer5[20],buffer6[20];

	for(int ii = 0; ii<m_nMs ; ii++) //Write x,y,vx,vy
	{
				

				//_gcvt( m_aM[ii].m_fX, 18, buffer1);
				//m_file.WriteString(buffer1);
				temp.Format("%lf",m_aM[ii].m_fX);
				m_file.WriteString(temp); 
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fY, 18, buffer2);
				//m_file.WriteString(buffer2); 
				temp.Format("%lf",m_aM[ii].m_fY);
				m_file.WriteString(temp); 
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fZ, 18, buffer5);
				//m_file.WriteString(buffer5); 
				temp.Format("%lf",m_aM[ii].m_fZ);
				m_file.WriteString(temp); 
				m_file.WriteString("\n"); 


				//_gcvt( m_aM[ii].m_fVx, 18, buffer3);
				//m_file.WriteString(buffer3);
				temp.Format("%lf",m_aM[ii].m_fVx);
				m_file.WriteString(temp); 
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fVy, 18, buffer4);
				//m_file.WriteString(buffer4); 
				temp.Format("%lf",m_aM[ii].m_fVy);
				m_file.WriteString(temp); 
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fVz, 18, buffer6);
				//m_file.WriteString(buffer6); 
				temp.Format("%lf",m_aM[ii].m_fVz);
				m_file.WriteString(temp); 
				m_file.WriteString("\n"); 

				m_file.WriteString("\n"); 
	}


	m_file.Close();

		m_file.Open("Heatflux.txt",CFile::modeCreate | CFile::modeWrite);

		//char heatflux1[20],heatflux2[20],heatflux3[20];

		//_gcvt( m_heatfluxx, 18, heatflux1);
		//m_file.WriteString(heatflux1);
		temp.Format("%.15e",m_heatfluxx);	m_file.WriteString(temp);
		m_file.WriteString("\n"); 

		//_gcvt( m_heatfluxy, 18, heatflux2);
		//m_file.WriteString(heatflux2);
		temp.Format("%.15e",m_heatfluxy);	m_file.WriteString(temp);
		m_file.WriteString("\n"); 

		//_gcvt( m_heatfluxz, 18, heatflux3);
		//m_file.WriteString(heatflux3);
		temp.Format("%.15e",m_heatfluxz);	m_file.WriteString(temp);
		m_file.WriteString("\n"); 



	m_file.Close();
	

}

void CMFCGLView::AssignImage(int ai, int bi, int ci, int ao, int bo, int co, double xx, double yy, double zz)
{
		for (int nn = 0; nn< Cell[ao][bo][co].NM ; nn++) 
		{	
		Cell[ai][bi][ci].x[nn] = Cell[ao][bo][co].x[nn] + xx;
		Cell[ai][bi][ci].y[nn] = Cell[ao][bo][co].y[nn] + yy;
		Cell[ai][bi][ci].z[nn] = Cell[ao][bo][co].z[nn] + zz;
		Cell[ai][bi][ci].mtype[nn] = 3; //set as image molecule
		Cell[ai][bi][ci].NM = Cell[ao][bo][co].NM;
		
		}
}

void CMFCGLView::AssignPerBC()
{

	int xc,yc,zc;
	double dNCX,dNCY,dNCZ;
	dNCX = m_fW/m_cellsize;
	dNCY = m_fH/m_cellsize;
	dNCZ = m_fZ/m_cellsize;

	// 6 face

	for (xc = 1; xc <= dNCX ; xc++)  //z = 0 plane
	{
		for ( yc = 1; yc <= dNCY ; yc++)
		{
			AssignImage (xc,yc,0,xc,yc,int(dNCZ),0,0,-m_fZ);
			AssignImage (xc,yc,int(dNCZ)+1,xc,yc,1,0,0,m_fZ);
		}
	}


	for (yc = 1; yc <= dNCY ; yc++)  //x = 0 plane
	{
		for ( zc = 1; zc <= dNCZ ; zc++)
		{
			AssignImage (0,yc,zc,int(dNCX),yc,zc,-m_fW,0,0);
			AssignImage (int(dNCX)+1,yc,zc,1,yc,zc,m_fW,0,0);
		}
	}

	//skip for wall case
	if (m_nwMs == 0)
	{
	for (zc = 1; zc <= dNCZ ; zc++)  //y = 0 plane
	{
		for ( xc = 1; xc <= dNCX ; xc++)
		{
			AssignImage (xc,0,zc,xc,int(dNCY),zc,0,-m_fH,0);
			AssignImage (xc,int(dNCY)+1,zc,xc,1,zc,0,m_fH,0);
		}
	}
	}


	//12 edge line
	
	//skip for wall case
	if (m_nwMs == 0)
	{
		//Rear-Top
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,int(dNCY+1),0,xc,1,int(dNCZ),0,m_fH,-m_fZ);}
		//Rear-Bottom
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,0,0,xc,int(dNCY),int(dNCZ),0,-m_fH,-m_fZ);}
		//Front-Bottom
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,0,int(dNCZ+1),xc,int(dNCY),1,0,-m_fH,m_fZ);}
		//Front-Top
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,int(dNCY+1),int(dNCZ+1),xc,1,1,0,m_fH,m_fZ);}
	}

	//Front-Left
	for (yc = 1; yc <= dNCY ; yc++) 
	{		AssignImage (0,yc,int(dNCZ+1),int(dNCX),yc,1,-m_fW,0,m_fZ);}
	//Front-Right
	for (yc = 1; yc <= dNCY ; yc++) 
	{		AssignImage (int(dNCX+1),yc,int(dNCZ+1),1,yc,1,m_fW,0,m_fZ);}
	//Rear-Left
	for (yc = 1; yc <= dNCY ; yc++) 
	{		AssignImage (0,yc,0,int(dNCX),yc,int(dNCZ),-m_fW,0,-m_fZ);}
	//Rear-Right
	for (yc = 1; yc <= dNCY ; yc++) 
	{		AssignImage (int(dNCX+1),yc,0,1,yc,int(dNCZ),m_fW,0,-m_fZ);}

		//skip for wall case
	if (m_nwMs == 0)
	{
		//Left-Bottom
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (0,0,zc,int(dNCX),int(dNCY),zc,-m_fW,-m_fH,0);}
		//Left-Top
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (0,int(dNCY+1),zc,int(dNCX),1,zc,-m_fW,m_fH,0);}
		//Right-Bottom
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (int(dNCX+1),0,zc,1,int(dNCY),zc,m_fW,-m_fH,0);}
		//Right-Top
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (int(dNCX+1),int(dNCY+1),zc,1,1,zc,m_fW,m_fH,0);}
	}

	//8 corner
		//skip for wall case
	if (m_nwMs == 0)
	{
		//Rear-Bottom-Left
		AssignImage (0,0,0,int(dNCX),int(dNCY),int(dNCZ),-m_fW,-m_fH,-m_fZ);
		//Rear-Bottom-Right
		AssignImage (int(dNCX)+1,0,0,1,int(dNCY),int(dNCZ),m_fW,-m_fH,-m_fZ);
		//Rear-Top-Left
		AssignImage (0,int(dNCY)+1,0,int(dNCX),1,int(dNCZ),-m_fW,m_fH,-m_fZ);
		//Rear-Top-Right
		AssignImage (int(dNCX)+1,int(dNCY)+1,0,1,1,int(dNCZ),m_fW,m_fH,-m_fZ);
		//============================================================================
		//Front-Top-Right
		AssignImage (int(dNCX)+1,int(dNCY)+1,int(dNCZ)+1,1,1,1,m_fW,m_fH,m_fZ);
		//Front-Top-Left
		AssignImage (0,int(dNCY)+1,int(dNCZ)+1,int(dNCX),1,1,-m_fW,m_fH,m_fZ);
		//Front-Bottom-Right
		AssignImage (int(dNCX)+1,0,int(dNCZ)+1,1,int(dNCY),1,m_fW,-m_fH,m_fZ);
		//Front-Bottom-Left
		AssignImage (0,0,int(dNCZ)+1,int(dNCX),int(dNCY),1,-m_fW,-m_fH,m_fZ);
	}

} //void ::AssignPerBC()


void CMFCGLView::SortWallMolecules()
{


		int ii;
		//m_dcutdist
		////For periodic boundary condition and sort for cell
		double dNCX,dNCY,dNCZ;
		dNCX = m_fW/m_cellsize;
		dNCY = m_fH/m_cellsize;
		dNCZ = m_fZ/m_cellsize;

		for( ii = 0; ii<m_nwMs ; ii++) // sort wall molecules
		{	


		//For X - periodic condition



			if (m_aTwM[ii].m_fXo >= m_fW) {	m_aTwM[ii].m_fXo = m_aTwM[ii].m_fXo - m_fW;
			m_aTwM[ii].m_fX = m_aTwM[ii].m_fX - m_fW;}
			if (m_aTwM[ii].m_fXo < 0.0)   {	m_aTwM[ii].m_fXo = m_aTwM[ii].m_fXo + m_fW;
			m_aTwM[ii].m_fX = m_aTwM[ii].m_fX + m_fW;	}

			if (m_aBwM[ii].m_fXo >= m_fW) {	m_aBwM[ii].m_fXo = m_aBwM[ii].m_fXo - m_fW;
			m_aBwM[ii].m_fX = m_aBwM[ii].m_fX - m_fW;}
			if (m_aBwM[ii].m_fXo < 0.0)   {	m_aBwM[ii].m_fXo = m_aBwM[ii].m_fXo + m_fW;	
			m_aBwM[ii].m_fX = m_aBwM[ii].m_fX + m_fW;	}


			//Sort into cells   			
			int xc,yc,zc;
            // zero cell is area for image
			//use original for cell put instant for position
			xc = int(m_aBwM[ii].m_fXo/m_cellsize)+1;  //decimal point removed.
			yc = int(m_aBwM[ii].m_fYo/m_cellsize);  //decimal point removed.
			zc = int(m_aBwM[ii].m_fZo/m_cellsize)+1;  //decimal point removed.
			
		//	if (xc < 0) xc =0; //if (xc >= int(dNCX)) xc = int(dNCX);
		//	if (yc < 0) yc =0; //if (yc >= int(dNCY)) yc = int(dNCY);
		//	if (zc < 0) zc =0; //if (zc >= int(dNCZ)) zc = int(dNCZ);

			Cell[xc][yc][zc].N[Cell[xc][yc][zc].NM] = ii; //Number of molecule - need to redistribute
			Cell[xc][yc][zc].x[Cell[xc][yc][zc].NM] = m_aBwM[ii].m_fX;
			Cell[xc][yc][zc].y[Cell[xc][yc][zc].NM] = m_aBwM[ii].m_fY;
			Cell[xc][yc][zc].z[Cell[xc][yc][zc].NM] = m_aBwM[ii].m_fZ;
			Cell[xc][yc][zc].mtype[Cell[xc][yc][zc].NM] = 2;
			if (m_aBwM[ii].m_fVy != 0.0) Cell[xc][yc][zc].mtype[Cell[xc][yc][zc].NM] = 4;
			Cell[xc][yc][zc].NM = Cell[xc][yc][zc].NM +1;
			

			

//		if (xc > int(dNCX) || zc > int(dNCZ) ) Cell[xc][yc][zc].NM = Cell[xc][yc][zc].NM -1;;
		//if (xc > int(dNCX) && zc > int(dNCZ) ) Cell[xc][yc][zc].NM = Cell[xc][yc][zc].NM +1;;
		

//use original for cell put instant for position
			xc = int(m_aTwM[ii].m_fXo/m_cellsize)+1;  //decimal point removed.
			yc = int(m_aTwM[ii].m_fYo/m_cellsize)+1;//+1;  //decimal point removed.
			zc = int(m_aTwM[ii].m_fZo/m_cellsize)+1;  //decimal point removed.

			if (yc >= int(dNCY)+1) yc = int(dNCY)+1;

		//	if (xc < 0) xc =0; //if (xc >= int(dNCX)) xc = int(dNCX);
		//	if (yc < 0) yc =0; 
		//	if (zc < 0) zc =0; //if (zc >= int(dNCZ)) zc = int(dNCZ);


			Cell[xc][yc][zc].N[Cell[xc][yc][zc].NM] = ii; //Number of molecule - need to redistribute
			Cell[xc][yc][zc].x[Cell[xc][yc][zc].NM] = m_aTwM[ii].m_fX;
			Cell[xc][yc][zc].y[Cell[xc][yc][zc].NM] = m_aTwM[ii].m_fY;
			Cell[xc][yc][zc].z[Cell[xc][yc][zc].NM] = m_aTwM[ii].m_fZ;
			Cell[xc][yc][zc].mtype[Cell[xc][yc][zc].NM] = 2;
			if (m_aTwM[ii].m_fVy != 0.0) Cell[xc][yc][zc].mtype[Cell[xc][yc][zc].NM] = 4;
			Cell[xc][yc][zc].NM = Cell[xc][yc][zc].NM +1;


			
//		if (xc > int(dNCX)|| zc > int(dNCZ) ) Cell[xc][yc][zc].NM = Cell[xc][yc][zc].NM -1;;
	
	  } 		//for(ii = 0; ii<m_nMs ; ii++)


	//duplication check
/*	
	for ( ii = 0 ; ii < m_nwMs  ; ii++)
		{
			for (int jj = 0 ; jj < m_nwMs  ; jj++)
			{
				if (ii!=jj)
				{
					if (m_aBwM[ii].m_fX == m_aBwM[jj].m_fX)
					{
						if (m_aBwM[ii].m_fY == m_aBwM[jj].m_fY)
						{
							if (m_aBwM[ii].m_fZ == m_aBwM[jj].m_fZ)
							{
								AfxMessageBox("duplicate");
		
							}
						}
					}
				}
			}//for (int jj = 0 ; jj < tem_number  ; jj++)
		}

*/

	//put wall images 

	int xc,zc;
	//8 edge line need change


	//skip for wall case
	if (m_nwMs != 0)
	{
		//Rear-Top
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,int(dNCY)+1,0  ,xc,int(dNCY+1),int(dNCZ)  ,0,0,-m_fZ);}
		//Rear-Bottom
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,0,0, xc,0,int(dNCZ), 0,0,-m_fZ);}
		//Front-Bottom
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,0,int(dNCZ)+1, xc,0,1, 0,0,m_fZ);}
		//Front-Top
		for (xc = 1; xc <= dNCX ; xc++) 
		{		AssignImage (xc,int(dNCY)+1,int(dNCZ)+1, xc,int(dNCY)+1,1, 0,0,m_fZ);}
	}


		//skip for wall case
	if (m_nwMs != 0)
	{
		//Left-Bottom
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (0,0,zc, int(dNCX),0,zc, -m_fW,0,0);}
		//Left-Top
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (0,int(dNCY+1),zc,int(dNCX),int(dNCY+1),zc,-m_fW,0,0);}
		//Right-Bottom
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (int(dNCX)+1,0,zc ,1,0,zc ,m_fW,0,0);}
		//Right-Top
		for (zc = 1; zc <= dNCZ ; zc++) 
		{		AssignImage (int(dNCX)+1,int(dNCY+1),zc,1,int(dNCY+1),zc,m_fW,0,0);}
	}



	//8 corner
		//skip for wall case
	if (m_nwMs != 0)
	{
		//Rear-Bottom-Left
		AssignImage (0,0,0,int(dNCX),0,int(dNCZ),-m_fW,0,-m_fZ);
		//Rear-Bottom-Right
		AssignImage (int(dNCX)+1,0,0, 1,0,int(dNCZ), m_fW,0,-m_fZ);
		//Rear-Top-Left
		AssignImage (0,int(dNCY)+1,0 ,int(dNCX),int(dNCY)+1,int(dNCZ),-m_fW,0,-m_fZ);
		//Rear-Top-Right
		AssignImage (int(dNCX)+1,int(dNCY)+1,0, 1,int(dNCY)+1,int(dNCZ),m_fW,0,-m_fZ);
		//============================================================================

		//Front-Top-Right
		AssignImage (int(dNCX)+1,int(dNCY)+1,int(dNCZ)+1, 1,int(dNCY)+1,1,m_fW,0,m_fZ);
		//Front-Top-Left
		AssignImage (0,int(dNCY)+1,int(dNCZ)+1,  int(dNCX),int(dNCY)+1,1,  -m_fW,0,m_fZ);
		//Front-Bottom-Right
		AssignImage (int(dNCX)+1,0,int(dNCZ)+1   ,1,0,1,    m_fW,0.0,m_fZ);
		//Front-Bottom-Left
		AssignImage (0,0,int(dNCZ)+1,    int(dNCX),0,1,      -m_fW,0,m_fZ);
	}


}


void CMFCGLView::SortGridInit()
{

	int xx,yy,zz,ii;
	//m_dcutdist
	////For periodic boundary condition and sort for cell
		double dNCX,dNCY,dNCZ;
		dNCX = m_fW/m_cellsize;
		dNCY = m_fH/m_cellsize;
		dNCZ = m_fZ/m_cellsize;
	
		for (xx = 0; xx <= int(dNCX)+1 ;xx++) //Initialize cell data
		{	
			for (yy = 0; yy <= int(dNCY)+1 ;yy++)
			{
				for (zz = 0; zz <= int(dNCZ)+1 ;zz++)
				{
					Cell[xx][yy][zz].NM = 0;
				}
			}
		}


		for( ii = 0; ii<m_nMs ; ii++) // sort fluid molecules
		{	
	
			//For X - periodic condition
			if (m_aM[ii].m_fX >= m_fW) {	m_aM[ii].m_fX = m_aM[ii].m_fX - m_fW;	}
			if (m_aM[ii].m_fX < 0.0)   {	m_aM[ii].m_fX = m_aM[ii].m_fX + m_fW;	}
			//For Z - periodic condition
			if (m_aM[ii].m_fZ >= m_fZ) {	m_aM[ii].m_fZ = m_aM[ii].m_fZ - m_fZ;	}
			if (m_aM[ii].m_fZ < 0.0)   {	m_aM[ii].m_fZ = m_aM[ii].m_fZ + m_fZ;	}
			//For Y - periodic condition // for the case of the all periodic conditions
			if (m_aM[ii].m_fY >= m_fH) {	m_aM[ii].m_fY = m_aM[ii].m_fY - m_fH;	}
			if (m_aM[ii].m_fY < 0.0)   {	m_aM[ii].m_fY = m_aM[ii].m_fY + m_fH;	}
			
			//Sort into cells   			
			int xc,yc,zc;
			xc = int(m_aM[ii].m_fX/m_cellsize)+1;  //decimal point removed.
			yc = int(m_aM[ii].m_fY/m_cellsize)+1;  //decimal point removed.
			zc = int(m_aM[ii].m_fZ/m_cellsize)+1;  //decimal point removed.

			//if (xc < 0) xc =0;
			//if (xc >= int(dNCX)-1) xc = int(dNCX) - 1;
			//if (yc < 0) yc =0;
			//if (yc >= int(dNCY)-1) yc = int(dNCY) - 1;

			Cell[xc][yc][zc].N[Cell[xc][yc][zc].NM] = ii; //Number of molecule - need to redistribute
			Cell[xc][yc][zc].x[Cell[xc][yc][zc].NM] = m_aM[ii].m_fX;
			Cell[xc][yc][zc].y[Cell[xc][yc][zc].NM] = m_aM[ii].m_fY;
			Cell[xc][yc][zc].z[Cell[xc][yc][zc].NM] = m_aM[ii].m_fZ;
			Cell[xc][yc][zc].mtype[Cell[xc][yc][zc].NM] = 1;
			Cell[xc][yc][zc].NM = Cell[xc][yc][zc].NM +1;
			
			//if I assign the periodic particle here, i will have a problem later
		} 		//for(ii = 0; ii<m_nMs ; ii++)

		//Sort fluid molecules to periodic

	AssignPerBC();

	//if there is walls, erase bottom and top image
		if (m_nwMs != 0) SortWallMolecules();



}

void CMFCGLView::CollectData()
{


		int ii,jj;


		int count[10];
//		int countfor_04[5000];

		for ( jj = 0; jj < 10 ; jj++)
		{	
			strG[jj] = "0";	
			dG[jj] = 0.0; 
	
			count[jj] = 0;
		}


		double pBV;
		pBV=0.0;

		for( ii = 0; ii<m_nMs ; ii++)
		{	for (jj = 0; jj < 10; jj++)
			{	if ( m_aM[ii].m_fY >= double(jj)*m_fH/10.0 && m_aM[ii].m_fY < double(jj+1)*m_fH/10.0)
				{	
					//double pBV;
					//pBV=0.0;
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 
					dG[jj] = dG[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
					count[jj] = count[jj] +1;
				}
			}
		}



		for (jj = 0; jj < 10; jj++)
		{	
			if (count[jj] != 0)
			{
			dG[jj] = (1.0/3.0)*m_Kb*dG[jj]/double(count[jj]);
			denCount[jj] = double(count[jj]);
			}
		}
		//Collecting velocity data -END


			////////////////////////////////////////////////////////////////////
		// Get Temperature
		//m_Kb - Mass divided by Boltzmann constant
		m_V  = 0.0;
		m_T  = 0.0;
		m_Pot = 0.0;
		m_P = 0;
		//Get potentilal and total velocity
		for( ii = 0; ii<m_nMs ; ii++)
		{
			if (m_gamma2==1) 
			{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 

			m_V = m_V +  pow(m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
			//for exclude flow direction 2/2
			//m_V = m_V +  pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
			//m_gamma = 2.0*m_Umax/m_fH;

			//m_Pot = m_Pot + m_vec.fPotSign * sqrt(pow(m_aM[ii].m_fPotx/1.0,2.0) + pow(m_aM[ii].m_fPoty/1.0,2.0));
			m_Pot =m_Pot + m_aM[ii].m_fPot;
			//m_V = m_V +  m_aM[ii].m_fVx*m_aM[ii].m_fVx + m_aM[ii].m_fVy*m_aM[ii].m_fVy;
			//m_Pot = m_Pot + m_vec.fPotSign * sqrt(m_aM[ii].m_fPotx*m_aM[ii].m_fPotx + m_aM[ii].m_fPoty*m_aM[ii].m_fPoty);
		} 


		//m_T = (1.0/2.0)*m_Kb*m_V/double(m_nMs);//for exclude flow direction 3/3
		m_T = (1.0/3.0)*m_Kb*m_V/double(m_nMs);


		m_Pot = 0.5*m_Pot*m_Kb/double(m_nMs);  // m_kb = m/kb because e is e/m

	
		m_Kin = 0.5*m_Kb*m_V/double(m_nMs);

		m_Tot = m_Kin + m_Pot; //Kin and Pot is
		
	// Let's Do Pressure Calculation!!   !! Need to check and do the  the validation Not completed.
		double virial;
		virial = 0;
		for( ii = 0; ii<m_nMs ; ii++)
		{
			virial = virial +  0.5*m_aM[ii].m_fStress;
		} 
		
		m_P = (-virial*m_Kb*m_BoltzN + 2.0*m_BoltzN*m_T)/(3*m_fH*m_fW*m_fZ); //need to validate with ref.
		// was ((1.0e-3 * 1.0e21 * 1.0e6)
		//m_p= g nm-1 ps-2
		//pascal = N/m2 = J/m3 = kg m2 s-2 / m3   = kg m-1 s-2
		//convert gnmos to pascal
		//((m_to_nm * s_to_ps* s_to_ps)/kg_to_g) = 9+12+12-3 = 30 
		m_P = m_P*((m_to_nm * s_to_ps* s_to_ps)/kg_to_g); //gnmps to Pa
		m_P = m_P*(1.0e-3); // unit pa to Kpa 

		//!! Need to check and do the  the validation
//	return;

		//*m_Kb*m_BoltzN//   for the force, the acceleration need to be multiplied by *m_Kb*m_BoltzN since e->e/m

		//


		//Get Wall Temperature
		double wvel,wvelsum,wvelsum1,wvelsum2,wvelsum3;
		int num1,num2,num3;
		num1=0;
		num2=0;
		num3=0;
		wvelsum = 0.0;
		wvelsum1 = 0.0;
		wvelsum2 = 0.0;
		wvelsum3 = 0.0;
		double tempy;
		for(ii = 0; ii<m_nwMs ; ii++)
		{	
				wvel = pow(m_aTwM[ii].m_fVx,2.0)+pow(m_aTwM[ii].m_fVy,2.0)+pow(m_aTwM[ii].m_fVz,2.0);
				wvelsum = wvelsum + wvel;

			tempy = m_fH+(m_fMdist/2.0);
			if (m_aTwM[ii].m_fYo - tempy < 0.0)
			{
				wvel = pow(m_aTwM[ii].m_fVx,2.0)+pow(m_aTwM[ii].m_fVy,2.0)+pow(m_aTwM[ii].m_fVz,2.0);
				wvelsum1 = wvelsum1 + wvel;
				num1++;
			}
			
			
			if (fabs(m_aTwM[ii].m_fYo - tempy) < m_fMdist/3.0)
			{
				wvel = pow(m_aTwM[ii].m_fVx,2.0)+pow(m_aTwM[ii].m_fVy,2.0)+pow(m_aTwM[ii].m_fVz,2.0);
				wvelsum2 = wvelsum2 + wvel;
				num2++;
			}

			
			if (m_aTwM[ii].m_fYo - tempy > 0.0)
			{
				wvel = pow(m_aTwM[ii].m_fVx,2.0)+pow(m_aTwM[ii].m_fVy,2.0)+pow(m_aTwM[ii].m_fVz,2.0);
				wvelsum3 = wvelsum3 + wvel;
				num3++;
			}
			
		}
		d_Twalltemp1 = (0.1/0.3)*m_Kb*wvelsum1/double(num1);
		d_Twalltemp2 = (0.1/0.3)*m_Kb*wvelsum2/double(num2);
		d_Twalltemp3 = (0.1/0.3)*m_Kb*wvelsum3/double(num3);
		d_Twalltemp = (0.1/0.3)*m_Kb*wvelsum/double(m_nwMs);


		wvelsum = 0.0;
		wvelsum1 = 0.0;
		wvelsum2 = 0.0;
		wvelsum3 = 0.0;
		num1=0;
		num2=0;
		num3=0;

			for(ii = 0; ii<m_nwMs ; ii++)
		{	
				wvel = pow(m_aBwM[ii].m_fVx,2.0)+pow(m_aBwM[ii].m_fVy,2.0)+pow(m_aBwM[ii].m_fVz,2.0);
				wvelsum = wvelsum + wvel;

			if (m_aBwM[ii].m_fYo == 0.0)
			{
				wvel = pow(m_aBwM[ii].m_fVx,2.0)+pow(m_aBwM[ii].m_fVy,2.0)+pow(m_aBwM[ii].m_fVz,2.0);
				wvelsum1 = wvelsum1 + wvel;
				num1++;
			}
			if (m_aBwM[ii].m_fYo == -(m_fMdist/2.0))
			{
				wvel = pow(m_aBwM[ii].m_fVx,2.0)+pow(m_aBwM[ii].m_fVy,2.0)+pow(m_aBwM[ii].m_fVz,2.0);
				wvelsum2 = wvelsum2 + wvel;
				num2++;
			}
			if (m_aBwM[ii].m_fYo == -m_fMdist)
			{
				wvel = pow(m_aBwM[ii].m_fVx,2.0)+pow(m_aBwM[ii].m_fVy,2.0)+pow(m_aBwM[ii].m_fVz,2.0);
				wvelsum3 = wvelsum3 + wvel;
				num3++;
			}
			
		}
		d_Bwalltemp1 = (0.1/0.3)*m_Kb*wvelsum1/double(num1);
		d_Bwalltemp2 = (0.1/0.3)*m_Kb*wvelsum2/double(num2);
		d_Bwalltemp3 = (0.1/0.3)*m_Kb*wvelsum3/double(num3);
		d_Bwalltemp = (0.1/0.3)*m_Kb*wvelsum/double(m_nwMs);


		////////////////////////////////////////////////////////////////////

}
//Heat flux calculation


void CMFCGLView::CollectEnergy()
{

	int ii,jj;
	double dkin;

	ii = 0;
	dkin = 0.0;

	//m_Kb  Mass divided by Boltzmann constant [argon]       
    //Epsylon - gas molecules   [nm]2/[ps]2 //Set this to e/m
					double pBV;
					pBV=0.0;
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV

		for( ii = 0; ii<m_nMs ; ii++)
		{
			//dkin = pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
			dkin = pow(m_aM[ii].m_fVx,2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
			m_aM[ii].m_dh = (0.5*m_Kb*dkin + 0.5*m_aM[ii].m_fPot*m_Kb)*m_BoltzN;	
			//unit for tis eqn = [g][nm]2/[ps]2  -> dimension for h same Dimension as [J]

			//pot times m_kb gives e/kb because epsylon is e/m
			//energy is J = N*m = kg*m^2/s^2  heat flux = W/m2   
			//J= kg m2 s-2
			//W= J s-1 = kg m2 s-3
			//[W] = [J]/[s] = [N][m]/[s] = [Kg]*[m]/[s^2]  * [m]/[s] = [Kg]*[m^2]/[s^3]
			// but here, 	m_aM[ii].m_dh is  g,nm,ps -> [g]*[nm^2]/[ps^2] energy
		}

		d_Jx = 0.0;
		d_Jy = 0.0;
		d_Jz = 0.0;

		for( ii = 0; ii<m_nMs ; ii++)
		{
			//d_Jx = d_Jx + m_aM[ii].m_dh*m_aM[ii].m_fVx - m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0);
			d_Jx = d_Jx + m_aM[ii].m_dh*m_aM[ii].m_fVx;
			d_Jy = d_Jy + m_aM[ii].m_dh*m_aM[ii].m_fVy;
			d_Jz = d_Jz + m_aM[ii].m_dh*m_aM[ii].m_fVz;

			//[g]*[nm^2]/[ps^2]  *  [nm]/[ps]

			//d_Jx = [g]*[nm^3]/[ps^3] = [watt] dimension * [Length] is what we are getting

			//[g]*[nm^2]/[ps^2] energy [J], d_J [g]*[nm^3]/[ps^3]
			
			//Thermal conductivity ramda or k? = [W m-2 K-1] 

			//Heat flux is defined as rate of heat transfer per unit cross-sectional area,
			//and is denoted q, resulting in units of watts per square metre

			//A heat current is a kinetic exchange rate between molecules,

			for( jj = 0; jj<m_nMs ; jj++)
			{
				if (ii != jj)
				{
				m_vec.fsignx = 0.0;	m_vec.fsigny = 0.0;	m_vec.accx = 0.0;	m_vec.accy = 0.0;
				m_vec.fForce = 0.0; m_vec.fPot = 0.0;
				GetReactionForce(m_aM[ii].m_fX,m_aM[jj].m_fX,m_aM[ii].m_fY,m_aM[jj].m_fY,m_aM[ii].m_fZ,m_aM[jj].m_fZ,m_dEgas);
				//m_vec.fsignx*m_vec.fForce*m_vec.accx   -> Forcex

				double forcex,forcey,forcez;
			
				forcex = m_vec.fsignx*m_vec.fForce*m_vec.accx*m_Kb*m_BoltzN;
				forcey = m_vec.fsigny*m_vec.fForce*m_vec.accy*m_Kb*m_BoltzN;
				forcez = m_vec.fsignz*m_vec.fForce*m_vec.accz*m_Kb*m_BoltzN;
				//m_aM[-].m_fax = m_aM[-].m_fax - m_vec.fsignx*m_vec.fForce*m_vec.accx/m_dMass_gas; //m_dMass_gas=1
				// m_kb = m/kb because epsylon is e/m

				//forcex = [nm]/[ps^2]    *   [g]
				//m_kb = m/kb because e is e/m   need to do something.
				//Epsylon - gas molecules   [nm]2/[ps]2 //Set this to e/m

				double dotFV;
				//dotFV = forcex*(m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)) + forcey*m_aM[ii].m_fVy + forcez*m_aM[ii].m_fVz;
				dotFV = forcex*m_aM[ii].m_fVx + forcey*m_aM[ii].m_fVy + forcez*m_aM[ii].m_fVz;

				d_Jx = d_Jx + 0.5*(m_aM[jj].m_fX - m_aM[ii].m_fX)*dotFV;
				d_Jy = d_Jy + 0.5*(m_aM[jj].m_fY - m_aM[ii].m_fY)*dotFV;
				d_Jz = d_Jz + 0.5*(m_aM[jj].m_fZ - m_aM[ii].m_fZ)*dotFV;
				//check unit [nm]  *  [kg]*[nm]/[ps^2]    *    [nm]/[ps]
				//=[kg]*[nm^3]/[ps^3]
				//d_Jx = [kg]*[nm^3]/[ps^3]
				//unit match!!!


				} //if (ii != jj)
			} //for( jj = 0; jj<m_nMs ; jj++) 
		}//for( ii = 0; ii<m_nMs ; ii++)

		//d_Jx = [g]*[nm^3]/[ps^3] = [watt] dimension * [Length] is what we are getting

		d_Jx = d_Jx/(m_fH*m_fW*m_fZ);
		d_Jy = d_Jy/(m_fH*m_fW*m_fZ);
		d_Jz = d_Jz/(m_fH*m_fW*m_fZ);

		//d_Jx = [g nm3 ps-3] [nm-3] = [g ps-3] = [g nm2 ps-3] [nm-2]
		//J= kg m2 s-2
		//W= J s-1 = kg m2 s-3, [W m-2] = [Kg s-3]

		//what we have is [g ps-3]
		//to change in watt/m2 (Kg s-3)

		d_Jx = d_Jx/(kg_to_g/(s_to_ps*s_to_ps*s_to_ps));
		d_Jy = d_Jy/(kg_to_g/(s_to_ps*s_to_ps*s_to_ps));
		d_Jz = d_Jz/(kg_to_g/(s_to_ps*s_to_ps*s_to_ps));

		//Now we have unit of watt / meter square!!!!!!

		//to get the thermal conductivity
///////////////////////  Get avaraged heat flux, ////////////////////////////////////////////

		if (d_ntimestep > 2500)
		{
		m_heatsumx = m_heatsumx + d_Jx;
		m_heatfluxx  = m_heatsumx/(d_ntimestep-2500.0);

		m_heatsumy = m_heatsumy + d_Jy;
		m_heatfluxy  = m_heatsumy/(d_ntimestep-2500.0);

		m_heatsumz = m_heatsumz + d_Jz;
		m_heatfluxz  = m_heatsumz/(d_ntimestep-2500.0);
		
		//m_heatfluxx  = d_Jx;  //heat
		//m_heatfluxy  = d_Jy;  //heat
		//m_heatfluxz  = d_Jz;  //heat
		}


///////////////////////////////Write output for time correlation function calc.		
		//Heat
/*		
		d_J[n_Jcount][0]=d_Jx;
		d_J[n_Jcount][1]=d_Jy;
		d_J[n_Jcount][2]=d_Jz;
		n_Jcount = n_Jcount + 1;


		if (n_Jcount == 100000)
		{
		char c_Jx[20],c_Jy[20],c_Jz[20],c_Jtotal[20];
		
			n_Jtotal = n_Jtotal + 1;
			_itoa(n_Jtotal,c_Jtotal,10);
			m_heatfluxfile.Open(c_Jtotal,CFile::modeCreate | CFile::modeWrite);
			//write all heat
				for (ii =0; ii<100000 ; ii++)
				{
				_gcvt(d_J[ii][0], 18, c_Jx);
				m_heatfluxfile.WriteString(c_Jx);
				m_heatfluxfile.WriteString("\n"); 
				_gcvt(d_J[ii][1], 18, c_Jy);
				m_heatfluxfile.WriteString(c_Jy);
				m_heatfluxfile.WriteString("\n"); 
				_gcvt(d_J[ii][2], 18, c_Jz);
				m_heatfluxfile.WriteString(c_Jz);
				m_heatfluxfile.WriteString("\n"); 
				d_J[ii][0]=0.0;
				d_J[ii][1]=0.0;
				d_J[ii][2]=0.0;
				}

			m_heatfluxfile.Close();
			n_Jcount = 0;
			TmpFileSave();

		}


//			m_ygridfile.Open("CollectV_5bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
//			char collectV_5bin_1st[20];
//			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
//			{
//					_gcvt(m_data_collectV5[ii], 18, collectV_5bin_1st);
//					m_ygridfile.WriteString(collectV_5bin_1st);
//					m_ygridfile.WriteString("\n"); 
//			}
//			m_ygridfile.Close();
*/
}


void CMFCGLView::DriveFlow()
{
			//For pressure driven case
		if (m_nSp == 1)
		{
			for(int ii = 0; ii<m_nMs ; ii++)
			{
				m_aM[ii].m_fax = m_aM[ii].m_fax + m_dGravP;//0.25;
			}
		}

		//For shear driven case, move upper wall to right
		if (m_nSp == 2)
		{
			
			
			
			for(int ii = 0; ii<m_nwMs ; ii++)
			{
			m_aTwM[ii].m_fX =m_aTwM[ii].m_fX + m_Umax*m_fdt;// m_dmove*3;
			m_aTwM[ii].m_fXo =m_aTwM[ii].m_fXo + m_Umax*m_fdt;// m_dmove*3;

			//m_dEgas = 4.0*0.02466367713004;    
			//Epsylon - gas molecules   [nm]2/[ps]2 //Set this to e/m
			//m_Umax = 0.15704673549628 nm/ps
				
			m_aBwM[ii].m_fX =m_aBwM[ii].m_fX - m_Umax*m_fdt;// m_dmove*3;
			m_aBwM[ii].m_fXo =m_aBwM[ii].m_fXo - m_Umax*m_fdt;// - m_dmove*3;
			}
		}

		//Vibrate wall
/*
		if (d_vibw != 0.0) d_vibw = d_vibw + 1.0;
		if (m_nWallDown == 1) 
		{
			for(int ii = 0; ii<m_nwMs ; ii++)
			{	
				//m_aTwM[ii].m_fY =m_aTwM[ii].m_fY + m_dmove*1.5;
				//m_aBwM[ii].m_fY =m_aBwM[ii].m_fY + m_dmove*3.0; //Wall Dow


				if  (d_vibw == 2)
					{
					double V1,V2,SD;
					double A1,A2;
					A1 = 0.0; A2 = 0.0;
					SD = sqrt((1.0/m_Kb)*m_Ttarget)*0.4;
					for (int rr=0;rr<12;rr++)
						{
							A1 = A1 + double(rand())/double(RAND_MAX);
							A2 = A2 + double(rand())/double(RAND_MAX);
						}
					V1 = ( A1 - 6.0 ) * SD;
					V2 = ( A2 - 6.0 ) * SD;
					m_V1[ii] = V1; 
					m_V2[ii] = V2;
					}
	
				//Wall Vibration
/*
				m_aTwM_old[ii].m_fX =m_aTwM[ii].m_fX;
				m_aBwM_old[ii].m_fX =m_aBwM[ii].m_fX;
				m_aTwM_old[ii].m_fY =m_aTwM[ii].m_fY;
				m_aBwM_old[ii].m_fY =m_aBwM[ii].m_fY;


				m_aTwM[ii].m_fX =m_aTwM_orig[ii].m_fX + 0.1*m_V1[ii]*sin((d_vibw+double(ii)-1.0)*0.2);
				m_aBwM[ii].m_fX =m_aBwM_orig[ii].m_fX + 0.1*m_V1[ii]*sin((d_vibw+double(ii)-1.0)*0.2); //Wall Down
				m_aTwM[ii].m_fY =m_aTwM_orig[ii].m_fY + 0.1*m_V2[ii]*sin((d_vibw+double(ii)-1.0)*0.2);
				m_aBwM[ii].m_fY =m_aBwM_orig[ii].m_fY + 0.1*m_V2[ii]*sin((d_vibw+double(ii)-1.0)*0.2); //Wall Down

				if (d_walltemp != 0.0)
				{
					double tempdiff;

				tempdiff = d_walltemp - 120.0;
				}

  
	
			}
			
		}
*/
 }

void CMFCGLView::GetCellInteraction()
{

	int xx,yy,zz,ii,jj,kk,nn,mm;
	double dNCX,dNCY,dNCZ;
	dNCX = m_fW/m_cellsize;
	dNCY = m_fH/m_cellsize;
	dNCZ = m_fZ/m_cellsize;

	//xx number of cell in x direction
	//yy number of cell in y direction
	//zz number of cell in z direction

	//Cell[xx][yy][zz].NM  Number of molecules in the cell
	//Cell[xx][yy][zz].N[ii] ii the molecule in the (xx,yy) cell

	// GetReactionForce(xi,xj,yi,yj,m_dEgas) 
	// m_dEgas = Epsylon / m
	// return a = a - m_vec.fsignx*m_vec.fForce*m_vec.accx;
int wallkey;
	wallkey=0;

//cell interaction for inner fluid
		
	
	for (xx = 1; xx < dNCX + 1 ; xx++)
	{
		for ( yy = 1; yy < dNCY + 1 ; yy++)
		{
			for ( zz = 1; zz < dNCZ + 1 ; zz++)
			{


						for (ii = -1; ii < 2 ; ii++)
						{
							for ( jj = -1; jj < 2 ; jj++)
							{
								for ( kk = -1; kk < 2 ; kk++)
								{

                                   // #pragma omp parallel for //reduction (+ : ax,ay,az,Pot,Stress,Virialxy)//still need some treatment for OpenMP //reduction 을 이용하면 많은수의 gasflow는 내부루프만 병렬처리하고 누적값 단일변수로 reduction 처리하면 가능할것같기도.
									for (nn = 0; nn< Cell[xx][yy][zz].NM ; nn++) 
									{
											for(mm = 0; mm < Cell[xx+ii][yy+jj][zz+kk].NM ; mm++) //left-middle
											{
											if ( ii==0 && jj == 0 && kk ==0 && nn==mm) continue;
												if (yy == 1 || yy == dNCY)
												{
													
													if (yy == 1 && jj == -1)
													{GetReactionForce(Cell[xx][yy][zz].x[nn],Cell[xx+ii][yy+jj][zz+kk].x[mm],Cell[xx][yy][zz].y[nn],Cell[xx+ii][yy+jj][zz+kk].y[mm],Cell[xx][yy][zz].z[nn],Cell[xx+ii][yy+jj][zz+kk].z[mm],m_dEwall);
													wallkey=1;}
													if (yy == dNCY && jj == 1)
													{GetReactionForce(Cell[xx][yy][zz].x[nn],Cell[xx+ii][yy+jj][zz+kk].x[mm],Cell[xx][yy][zz].y[nn],Cell[xx+ii][yy+jj][zz+kk].y[mm],Cell[xx][yy][zz].z[nn],Cell[xx+ii][yy+jj][zz+kk].z[mm],m_dEwall);
													wallkey=1;}
													if (wallkey == 0)
													{
														GetReactionForce(Cell[xx][yy][zz].x[nn],Cell[xx+ii][yy+jj][zz+kk].x[mm],Cell[xx][yy][zz].y[nn],Cell[xx+ii][yy+jj][zz+kk].y[mm],Cell[xx][yy][zz].z[nn],Cell[xx+ii][yy+jj][zz+kk].z[mm],m_dEgas);
													}
													wallkey=0;
												}
												else
												{
												GetReactionForce(Cell[xx][yy][zz].x[nn],Cell[xx+ii][yy+jj][zz+kk].x[mm],Cell[xx][yy][zz].y[nn],Cell[xx+ii][yy+jj][zz+kk].y[mm],Cell[xx][yy][zz].z[nn],Cell[xx+ii][yy+jj][zz+kk].z[mm],m_dEgas);
												}
											
											//need to change here for OPENMP static data member is acceptable
											m_aM[Cell[xx][yy][zz].N[nn]].m_fax = m_aM[Cell[xx][yy][zz].N[nn]].m_fax + m_vec.fsignx*m_vec.fForce*m_vec.accx; //m_dMass_gas=1
											m_aM[Cell[xx][yy][zz].N[nn]].m_fay = m_aM[Cell[xx][yy][zz].N[nn]].m_fay + m_vec.fsigny*m_vec.fForce*m_vec.accy;
											m_aM[Cell[xx][yy][zz].N[nn]].m_faz = m_aM[Cell[xx][yy][zz].N[nn]].m_faz + m_vec.fsignz*m_vec.fForce*m_vec.accz;
											m_aM[Cell[xx][yy][zz].N[nn]].m_fPot =m_aM[Cell[xx][yy][zz].N[nn]].m_fPot + m_vec.fPot;
											m_aM[Cell[xx][yy][zz].N[nn]].m_fStress = m_aM[Cell[xx][yy][zz].N[nn]].m_fStress + m_vec.fdist*m_vec.fForce;
											//0.5*rxij*Fyij
											m_aM[Cell[xx][yy][zz].N[nn]].m_fVirialxy = m_aM[Cell[xx][yy][zz].N[nn]].m_fVirialxy + (m_vec.fdist*m_vec.fsignx*m_vec.accx)*(m_vec.fsigny*m_vec.fForce*m_vec.accy);

											//#pragma omp atomic //Reserve for OpenMP
											//m_aM_OMP[Cell[xx][yy][zz].N[nn]].m_fax += m_vec.fsignx*m_vec.fForce*m_vec.accx; //m_dMass_gas=1
											//m_aM_OMP[Cell[xx][yy][zz].N[nn]].m_fay += m_vec.fsigny*m_vec.fForce*m_vec.accy;
											//m_aM_OMP[Cell[xx][yy][zz].N[nn]].m_faz += m_vec.fsignz*m_vec.fForce*m_vec.accz;
											//m_aM_OMP[Cell[xx][yy][zz].N[nn]].m_fPot += m_vec.fPot;
											//m_aM_OMP[Cell[xx][yy][zz].N[nn]].m_fStress += m_vec.fdist*m_vec.fForce;//	0.5*rxij*Fyij
											//m_aM_OMP[Cell[xx][yy][zz].N[nn]].m_fVirialxy += (m_vec.fdist*m_vec.fsignx*m_vec.accx)*(m_vec.fsigny*m_vec.fForce*m_vec.accy);
											}

									}
									



								}
							}
						}


			}
		}
	}



// for wall interaction
	//!need to add - if wall temp != 0 
	if (m_nwMs != 0)
	{
		
			yy=0;  //for bottom wall
			for (xx = 1; xx < dNCX + 1 ; xx++)  // 1 cell to DNCX cell  (inner fluid cell) 
			{
					for ( zz = 1; zz < dNCZ + 1 ; zz++)
					{

								jj=1;
								for (ii = -1; ii < 2 ; ii++)   //-1 0 1  (for surroundiung calc)
								{
										for ( kk = -1; kk < 2 ; kk++)
										{
											for (nn = 0; nn< Cell[xx][yy][zz].NM ; nn++) 
											{
												for(mm = 0; mm < Cell[xx+ii][yy+jj][zz+kk].NM ; mm++) //left-middle
												{
												//	if ( ii==0 && jj == 0 && kk ==0 && nn==mm) continue;  //skip for the same molecules
												GetReactionForce(Cell[xx][yy][zz].x[nn],Cell[xx+ii][yy+jj][zz+kk].x[mm],Cell[xx][yy][zz].y[nn],Cell[xx+ii][yy+jj][zz+kk].y[mm],Cell[xx][yy][zz].z[nn],Cell[xx+ii][yy+jj][zz+kk].z[mm],m_dEwall);
												m_aBwM[Cell[xx][yy][zz].N[nn]].m_fax = m_aBwM[Cell[xx][yy][zz].N[nn]].m_fax + m_vec.fsignx*m_vec.fForce*m_vec.accx;
												m_aBwM[Cell[xx][yy][zz].N[nn]].m_fay = m_aBwM[Cell[xx][yy][zz].N[nn]].m_fay + m_vec.fsigny*m_vec.fForce*m_vec.accy;
												m_aBwM[Cell[xx][yy][zz].N[nn]].m_faz = m_aBwM[Cell[xx][yy][zz].N[nn]].m_faz + m_vec.fsignz*m_vec.fForce*m_vec.accz;
												m_aBwM[Cell[xx][yy][zz].N[nn]].m_fPot =m_aBwM[Cell[xx][yy][zz].N[nn]].m_fPot + m_vec.fPot;
												m_aBwM[Cell[xx][yy][zz].N[nn]].m_fStress = m_aBwM[Cell[xx][yy][zz].N[nn]].m_fStress + m_vec.fdist*m_vec.fForce;
												}
											}
										}//for ( kk = -1; kk < 2 ; kk++)
								}//for (ii = -1; ii < 2 ; ii++) 


					}//	for ( zz = 1; zz < dNCZ + 1 ; zz++)
			}//for (xx = 1; xx < dNCX + 1 ; xx++)  // 1 cell to DNCX cell  (inner fluid cell) 


			yy= int(dNCY)+1;  //for bottom wall
			for (xx = 1; xx < dNCX + 1 ; xx++)  // 1 cell to DNCX cell  (inner fluid cell) 
			{
					for ( zz = 1; zz < dNCZ + 1 ; zz++)
					{

								jj=-1;
								for (ii = -1; ii < 2 ; ii++)   //-1 0 1  (for surroundiung calc)
								{
										for ( kk = -1; kk < 2 ; kk++)
										{
											for (nn = 0; nn< Cell[xx][yy][zz].NM ; nn++) 
											{
												for(mm = 0; mm < Cell[xx+ii][yy+jj][zz+kk].NM ; mm++) //left-middle
												{
												//	if ( ii==0 && jj == 0 && kk ==0 && nn==mm) continue;  //skip for the same molecules
												GetReactionForce(Cell[xx][yy][zz].x[nn],Cell[xx+ii][yy+jj][zz+kk].x[mm],Cell[xx][yy][zz].y[nn],Cell[xx+ii][yy+jj][zz+kk].y[mm],Cell[xx][yy][zz].z[nn],Cell[xx+ii][yy+jj][zz+kk].z[mm],m_dEwall);
												m_aTwM[Cell[xx][yy][zz].N[nn]].m_fax = m_aTwM[Cell[xx][yy][zz].N[nn]].m_fax + m_vec.fsignx*m_vec.fForce*m_vec.accx;
												m_aTwM[Cell[xx][yy][zz].N[nn]].m_fay = m_aTwM[Cell[xx][yy][zz].N[nn]].m_fay + m_vec.fsigny*m_vec.fForce*m_vec.accy;
												m_aTwM[Cell[xx][yy][zz].N[nn]].m_faz = m_aTwM[Cell[xx][yy][zz].N[nn]].m_faz + m_vec.fsignz*m_vec.fForce*m_vec.accz;
												m_aTwM[Cell[xx][yy][zz].N[nn]].m_fPot =m_aTwM[Cell[xx][yy][zz].N[nn]].m_fPot + m_vec.fPot;
												m_aTwM[Cell[xx][yy][zz].N[nn]].m_fStress = m_aTwM[Cell[xx][yy][zz].N[nn]].m_fStress + m_vec.fdist*m_vec.fForce;
												}
											}
										}//for ( kk = -1; kk < 2 ; kk++)
								}//for (ii = -1; ii < 2 ; ii++) 


					}//	for ( zz = 1; zz < dNCZ + 1 ; zz++)
			}//for (xx = 1; xx < dNCX + 1 ; xx++)  // 1 cell to DNCX cell  (inner fluid cell) 




	}//	if (m_nwMs != 0)

	//	m_aTwM[ii].m_fX
	//double m_fX,m_fY,m_fZ;    //Position
	//double m_fXo,m_fYo,m_fZo; //Original Position



//need to implement bottom wall and top wall


} //-End of GetCellInteraction() 


void CMFCGLView::GetReactionForce(double xi,double xj, double yi,double yj, double zi,double zj,double e)
{


		//double fang;
		//fang =double(m_nMs) 0.0;

		m_vec.fsignx = 0.0;	m_vec.fsigny = 0.0;	m_vec.accx = 0.0; m_vec.accy = 0.0;
		m_vec.fForce = 0.0;	m_vec.fPot = 0.0; m_vec.fdist = 0.0;


		m_vec.fdist = sqrt(pow((xi - xj),2.0)+pow((yi - yj),2.0)+pow((zi - zj),2.0));

		if (m_vec.fdist > m_dcutdist) return; //cut off

		//if (m_vec.fdist < 0.5)
		//{
		//	return; //cut off
		//}

		//m_vec.fdist = 2.5*m_dSigma;
        //m_vec.fForce = 0.0028248767111193 //need to fix
		double FSwitch;
		//FSwitch = 0.00036144369804767;
		 FSwitch = 0.00053721927743507;
		//FSwitch = 0;
		// cutdist = 0.382*3.0; //= 1.143   // FSwitch 0.00036144369804767

		//m_vec.fdist = 1.143;
		
		//m_vec.fdist = 0.381;  //-0.015417612670157
	//	m_vec.fdist = 0.382;  //-0.0027856113856136   -0.0029613869650010
	//	m_vec.fdist = 0.383;  //0.0091721719477844   0.0089963963683970
		
//		m_vec.fdist = 0.540*2;//	0.00053721927743507

		m_vec.fForce = -e*(-12.0*(1.0/m_dSigma)*pow(m_dSigma/(m_vec.fdist),13.0)+6.0*(1.0/m_dSigma)*pow(m_dSigma/(m_vec.fdist),7.0)) + FSwitch;
		//force is -dV/vr 
     		 
				//get Potential 
		//need to compensate when pot comes into play. -> Tompson shear paper reference -> Clemson Univ book
		m_vec.fPot = e*(pow(m_dSigma/(m_vec.fdist),12.0)-pow(m_dSigma/(m_vec.fdist),6.0));
				
		if ( (xi - xj) > 0.0)	{m_vec.fsignx = 1.0;}
		if ( (xi - xj) < 0.0)	{m_vec.fsignx = -1.0;}
		if ( (xi - xj) == 0.0)	{m_vec.fsignx = 0.0;}

		if ( (yi - yj) > 0.0)	{m_vec.fsigny = 1.0;}
		if ( (yi - yj) < 0.0)	{m_vec.fsigny = -1.0;}
		if ( (yi - yj) == 0.0)	{m_vec.fsigny = 0.0;}

		if ( (zi - zj) > 0.0)	{m_vec.fsignz = 1.0;}
		if ( (zi - zj) < 0.0)	{m_vec.fsignz = -1.0;}
		if ( (zi - zj) == 0.0)	{m_vec.fsignz = 0.0;}
		
		//m_aM[].m_fay = m_aM[].m_fay - m_vec.fsigny*m_vec.fForce*m_vec.accy;
		
		//1.  m_fVirialxy = m_fVirialxy + 0.5*(m_vec.fsignx*m_vec.accx*m_vec.fdist)*(m_vec.fsigny*m_vec.accy*m_vec.fForce);
		//2.  m_fVirialxy = m_fVirialxy - 0.5*(m_vec.fsignx*m_vec.accx*m_vec.fdist)*(m_vec.fsigny*m_vec.accy*m_vec.fForce);


//				if ((xi - xj) != 0.0)
//				{
//					if ((yi - yj) != 0.0)
//					{
//						if ((zi - zj) != 0.0)
//						{
		m_vec.accx = (xi - xj)/m_vec.fdist;
		m_vec.accy = (yi - yj)/m_vec.fdist;
		m_vec.accz = (zi - zj)/m_vec.fdist;
		if (m_vec.accx <0.0) m_vec.accx = -m_vec.accx;			
		if (m_vec.accy <0.0) m_vec.accy = -m_vec.accy;
		if (m_vec.accz <0.0) m_vec.accz = -m_vec.accz;
//						}
//					}
//				}
/*
				//m_vec.accy = sin(fang);
				if ((yi - yj) == 0.0)
				{
				m_vec.accy = 0.0;// m_vec.accx = 1.0;
				}		
				
				//m_vec.accx = cos(fang);
				if ((xi - xj) == 0.0)
				{
				m_vec.accx = 0.0;// m_vec.accy = 1.0;
				}

				if ((zi - zj) == 0.0)
				{
				m_vec.accz = 0.0;// m_vec.accz = 1.0;
				}
*/



}

void CMFCGLView::CollectAndSaveFile()
{
	int ii,jj;
	

//Collect Data For long time averaging collection

		//Initialize
		for (ii =0; ii<500 ; ii++)
		{
			m_data_count5[ii] = 0.0;	 m_data_count10[ii] = 0.0;
			m_data_count50[ii] = 0.0;	 m_data_count100[ii] = 0.0;
			m_data_count500[ii] = 0.0;

			m_data_temp5[ii] = 0.0;	 m_data_temp10[ii] = 0.0;
			m_data_temp50[ii] = 0.0;	 m_data_temp100[ii] = 0.0;
			m_data_temp500[ii] = 0.0;

			m_data_vel5[ii] = 0.0;	 m_data_vel10[ii] = 0.0;
			m_data_vel50[ii] = 0.0;	 m_data_vel100[ii] = 0.0;
			m_data_vel500[ii] = 0.0;
		
			m_data_virial5[ii] = 0.0;	 m_data_virial10[ii] = 0.0;
			m_data_stress5[ii] = 0.0;	 m_data_stress10[ii] = 0.0;
			
			m_data_virial100[ii] = 0.0;
			m_data_stress100[ii] = 0.0;
		}


		double pBV;
		pBV=0.0;


		for( ii = 0; ii<m_nMs ; ii++)
		{

			for (jj = 0; jj < 5; jj++)
			{	if ( m_aM[ii].m_fY >= double(jj)*(m_fH/5.0) && m_aM[ii].m_fY < double(jj+1)*(m_fH/5.0))
				{	
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 

					m_data_count5[jj] = m_data_count5[jj] +1;  //density profile
					m_data_temp5[jj] = m_data_temp5[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
					//m_gamma = 2.0*m_Umax/m_fH;
					m_data_vel5[jj] = m_data_vel5[jj] + m_aM[ii].m_fVx;

					m_data_virial5[jj] = m_data_virial5[jj] + (m_aM[ii].m_fVx*m_aM[ii].m_fVy+ 0.5*m_aM[ii].m_fVirialxy)*m_Kb*m_BoltzN;
					m_data_stress5[jj] = m_data_stress5[jj] + (0.5*m_aM[ii].m_fVirialxy)*m_Kb*m_BoltzN;
				}
			}


			for (jj = 0; jj < 10; jj++)
			{	if ( m_aM[ii].m_fY >= double(jj)*(m_fH/10.0) && m_aM[ii].m_fY < double(jj+1)*(m_fH/10.0))
				{	
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 

					m_data_count10[jj] = m_data_count10[jj] +1;  //density profile
					m_data_temp10[jj] = m_data_temp10[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
					m_data_vel10[jj] = m_data_vel10[jj] + m_aM[ii].m_fVx;
					
					m_data_virial10[jj] = m_data_virial10[jj] + (m_aM[ii].m_fVx*m_aM[ii].m_fVy+ 0.5*m_aM[ii].m_fVirialxy)*m_Kb*m_BoltzN;
					m_data_stress10[jj] = m_data_stress10[jj] + (0.5*m_aM[ii].m_fVirialxy)*m_Kb*m_BoltzN;
				}
			}

			for (jj = 0; jj < 50; jj++)
			{	if ( m_aM[ii].m_fY >= double(jj)*(m_fH/40.0) && m_aM[ii].m_fY < double(jj+1)*(m_fH/40.0))
				{	
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 

					m_data_count50[jj] = m_data_count50[jj] +1;  //density profile
					m_data_temp50[jj] = m_data_temp50[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);;
					m_data_vel50[jj] = m_data_vel50[jj] + m_aM[ii].m_fVx;
				}
			}

			for (jj = 0; jj < 100; jj++)
			{	if ( m_aM[ii].m_fY >= double(jj)*(m_fH/100.0) && m_aM[ii].m_fY < double(jj+1)*(m_fH/100.0))
				{	
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 

					m_data_count100[jj] = m_data_count100[jj] +1;  //density profile
					m_data_temp100[jj] = m_data_temp100[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);;
					m_data_vel100[jj] = m_data_vel100[jj] + m_aM[ii].m_fVx;

					m_data_virial100[jj] = m_data_virial100[jj] + (m_aM[ii].m_fVx*m_aM[ii].m_fVy+ 0.5*m_aM[ii].m_fVirialxy)*m_Kb*m_BoltzN;
					m_data_stress100[jj] = m_data_stress100[jj] + (0.5*m_aM[ii].m_fVirialxy)*m_Kb*m_BoltzN;
				}
			}

			for (jj = 0; jj < 25; jj++)
			{	if ( m_aM[ii].m_fY >= double(jj)*(m_fH/25.0) && m_aM[ii].m_fY < double(jj+1)*(m_fH/25.0))
				{	
					if (m_gamma2==1) 
					{pBV = ((m_aM[ii].m_fY/m_fH)*(1.0-m_aM[ii].m_fY/m_fH)*0.0546*4.0);} //-m_gamma2*pBV 

					m_data_count500[jj] = m_data_count500[jj] +1;  //density profile 
					m_data_temp500[jj] = m_data_temp500[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0)  + pow(m_aM[ii].m_fVz,2.0);;
					m_data_vel500[jj] = m_data_vel500[jj] + m_aM[ii].m_fVx;
				}
			}


		}

		
			for (jj = 0; jj < 5; jj++)
			{				
				if (m_data_count5[jj] != 0.0)
				{m_data_temp5[jj] = (1.0/3.0)*m_Kb*m_data_temp5[jj]/m_data_count5[jj];
				m_data_vel5[jj] = m_data_vel5[jj]/m_data_count5[jj];}
				else{m_data_temp5[jj] = 0.0;
					m_data_vel5[jj]=0.0;}
			}

			for (jj = 0; jj < 10; jj++)
			{				
				if (m_data_count10[jj] != 0.0)
				{m_data_temp10[jj] = (1.0/3.0)*m_Kb*m_data_temp10[jj]/m_data_count10[jj];
				m_data_vel10[jj] = m_data_vel10[jj]/m_data_count10[jj];}
				else{m_data_temp10[jj] = 0.0;
					m_data_vel10[jj]=0.0;}
			}

			for (jj = 0; jj < 40; jj++)
			{				
				if (m_data_count50[jj] != 0.0)
				{m_data_temp50[jj] = (1.0/3.0)*m_Kb*m_data_temp50[jj]/m_data_count50[jj];
				m_data_vel50[jj] = m_data_vel50[jj]/m_data_count50[jj];}
				else{m_data_temp50[jj] = 0.0;
					m_data_vel50[jj]=0.0;}
			}
			for (jj = 0; jj < 100; jj++)
			{				
				if (m_data_count100[jj] != 0.0)
				{m_data_temp100[jj] = (1.0/3.0)*m_Kb*m_data_temp100[jj]/m_data_count100[jj];
				m_data_vel100[jj] = m_data_vel100[jj]/m_data_count100[jj];}
				else{m_data_temp100[jj] = 0.0;
				m_data_vel100[jj]=0.0;}
			}
			for (jj = 0; jj < 25; jj++)
			{				
				if (m_data_count500[jj] != 0.0)
				{m_data_temp500[jj] = (1.0/3.0)*m_Kb*m_data_temp500[jj]/m_data_count500[jj];
				m_data_vel500[jj] = m_data_vel500[jj]/m_data_count500[jj];}
				else{m_data_temp500[jj] = 0.0;
				m_data_vel500[jj]=0.0;}
			}


		if (d_ntimestep > m_startstep1 && d_ntimestep < m_endstep1+1)
		{ 
			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = m_data_collectD5[jj] + m_data_count5[jj]/double(m_averagingstep1);
			 m_data_collectT5[jj] = m_data_collectT5[jj] + m_data_temp5[jj]/double(m_averagingstep1);
			 m_data_collectV5[jj] = m_data_collectV5[jj] + m_data_vel5[jj]/double(m_averagingstep1);
			m_data_collectstress5[jj] = m_data_collectstress5[jj] + m_data_stress5[jj]/double(m_averagingstep1);
			m_data_collectvirial5[jj] = m_data_collectvirial5[jj] + m_data_virial5[jj]/double(m_averagingstep1);}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = m_data_collectD10[jj] + m_data_count10[jj]/double(m_averagingstep1);
			 m_data_collectT10[jj] = m_data_collectT10[jj] + m_data_temp10[jj]/double(m_averagingstep1);
			m_data_collectV10[jj] = m_data_collectV10[jj] + m_data_vel10[jj]/double(m_averagingstep1);
			m_data_collectstress10[jj] = m_data_collectstress10[jj] + m_data_stress10[jj]/double(m_averagingstep1);
			m_data_collectvirial10[jj] = m_data_collectvirial10[jj] + m_data_virial10[jj]/double(m_averagingstep1);}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = m_data_collectD50[jj] + m_data_count50[jj]/double(m_averagingstep1);
			 m_data_collectT50[jj] = m_data_collectT50[jj] + m_data_temp50[jj]/double(m_averagingstep1);
			 m_data_collectV50[jj] = m_data_collectV50[jj] + m_data_vel50[jj]/double(m_averagingstep1);}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = m_data_collectD100[jj] + m_data_count100[jj]/double(m_averagingstep1);
			 m_data_collectT100[jj] = m_data_collectT100[jj] + m_data_temp100[jj]/double(m_averagingstep1);
			m_data_collectV100[jj] = m_data_collectV100[jj] + m_data_vel100[jj]/double(m_averagingstep1);
			m_data_collectstress100[jj] = m_data_collectstress100[jj] + m_data_stress100[jj]/double(m_averagingstep1);
			m_data_collectvirial100[jj] = m_data_collectvirial100[jj] + m_data_virial100[jj]/double(m_averagingstep1);}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = m_data_collectD500[jj] + m_data_count500[jj]/double(m_averagingstep1);
			 m_data_collectT500[jj] = m_data_collectT500[jj] + m_data_temp500[jj]/double(m_averagingstep1);
			m_data_collectV500[jj] = m_data_collectV500[jj] + m_data_vel500[jj]/double(m_averagingstep1);}
		}

		if (d_ntimestep > m_startstep2 && d_ntimestep < m_endstep2+1)
		{ 
			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = m_data_collectD5[jj] + m_data_count5[jj]/double(m_averagingstep2);
			 m_data_collectT5[jj] = m_data_collectT5[jj] + m_data_temp5[jj]/double(m_averagingstep2);
			 m_data_collectV5[jj] = m_data_collectV5[jj] + m_data_vel5[jj]/double(m_averagingstep2);
			 m_data_collectstress5[jj] = m_data_collectstress5[jj] + m_data_stress5[jj]/double(m_averagingstep2);
			m_data_collectvirial5[jj] = m_data_collectvirial5[jj] + m_data_virial5[jj]/double(m_averagingstep2);}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = m_data_collectD10[jj] + m_data_count10[jj]/double(m_averagingstep2);
			 m_data_collectT10[jj] = m_data_collectT10[jj] + m_data_temp10[jj]/double(m_averagingstep2);
			m_data_collectV10[jj] = m_data_collectV10[jj] + m_data_vel10[jj]/double(m_averagingstep2);
			m_data_collectstress10[jj] = m_data_collectstress10[jj] + m_data_stress10[jj]/double(m_averagingstep2);
			m_data_collectvirial10[jj] = m_data_collectvirial10[jj] + m_data_virial10[jj]/double(m_averagingstep2);}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = m_data_collectD50[jj] + m_data_count50[jj]/double(m_averagingstep2);
			 m_data_collectT50[jj] = m_data_collectT50[jj] + m_data_temp50[jj]/double(m_averagingstep2);
			 m_data_collectV50[jj] = m_data_collectV50[jj] + m_data_vel50[jj]/double(m_averagingstep2);}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = m_data_collectD100[jj] + m_data_count100[jj]/double(m_averagingstep2);
			 m_data_collectT100[jj] = m_data_collectT100[jj] + m_data_temp100[jj]/double(m_averagingstep2);
			m_data_collectV100[jj] = m_data_collectV100[jj] + m_data_vel100[jj]/double(m_averagingstep2);
			m_data_collectstress100[jj] = m_data_collectstress100[jj] + m_data_stress100[jj]/double(m_averagingstep2);
			m_data_collectvirial100[jj] = m_data_collectvirial100[jj] + m_data_virial100[jj]/double(m_averagingstep2);}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = m_data_collectD500[jj] + m_data_count500[jj]/double(m_averagingstep2);
			 m_data_collectT500[jj] = m_data_collectT500[jj] + m_data_temp500[jj]/double(m_averagingstep2);
			m_data_collectV500[jj] = m_data_collectV500[jj] + m_data_vel500[jj]/double(m_averagingstep2);}
		}

		if (d_ntimestep > m_startstep3 && d_ntimestep < m_endstep3+1)
		{ 
			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = m_data_collectD5[jj] + m_data_count5[jj]/double(m_averagingstep3);
			 m_data_collectT5[jj] = m_data_collectT5[jj] + m_data_temp5[jj]/double(m_averagingstep3);
			 m_data_collectV5[jj] = m_data_collectV5[jj] + m_data_vel5[jj]/double(m_averagingstep3);
			 m_data_collectstress5[jj] = m_data_collectstress5[jj] + m_data_stress5[jj]/double(m_averagingstep3);
			m_data_collectvirial5[jj] = m_data_collectvirial5[jj] + m_data_virial5[jj]/double(m_averagingstep3);}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = m_data_collectD10[jj] + m_data_count10[jj]/double(m_averagingstep3);
			 m_data_collectT10[jj] = m_data_collectT10[jj] + m_data_temp10[jj]/double(m_averagingstep3);
			m_data_collectV10[jj] = m_data_collectV10[jj] + m_data_vel10[jj]/double(m_averagingstep3);
			m_data_collectstress10[jj] = m_data_collectstress10[jj] + m_data_stress10[jj]/double(m_averagingstep3);
			m_data_collectvirial10[jj] = m_data_collectvirial10[jj] + m_data_virial10[jj]/double(m_averagingstep3);}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = m_data_collectD50[jj] + m_data_count50[jj]/double(m_averagingstep3);
			 m_data_collectT50[jj] = m_data_collectT50[jj] + m_data_temp50[jj]/double(m_averagingstep3);
			 m_data_collectV50[jj] = m_data_collectV50[jj] + m_data_vel50[jj]/double(m_averagingstep3);}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = m_data_collectD100[jj] + m_data_count100[jj]/double(m_averagingstep3);
			 m_data_collectT100[jj] = m_data_collectT100[jj] + m_data_temp100[jj]/double(m_averagingstep3);
			m_data_collectV100[jj] = m_data_collectV100[jj] + m_data_vel100[jj]/double(m_averagingstep3);
			m_data_collectstress100[jj] = m_data_collectstress100[jj] + m_data_stress100[jj]/double(m_averagingstep3);
			m_data_collectvirial100[jj] = m_data_collectvirial100[jj] + m_data_virial100[jj]/double(m_averagingstep3);}

			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = m_data_collectD500[jj] + m_data_count500[jj]/double(m_averagingstep3);
			 m_data_collectT500[jj] = m_data_collectT500[jj] + m_data_temp500[jj]/double(m_averagingstep3);
			m_data_collectV500[jj] = m_data_collectV500[jj] + m_data_vel500[jj]/double(m_averagingstep3);}
		}

		if (d_ntimestep > m_startstep4 && d_ntimestep < m_endstep4+1)
		{ 
			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = m_data_collectD5[jj] + m_data_count5[jj]/double(m_averagingstep4);
			 m_data_collectT5[jj] = m_data_collectT5[jj] + m_data_temp5[jj]/double(m_averagingstep4);
			 m_data_collectV5[jj] = m_data_collectV5[jj] + m_data_vel5[jj]/double(m_averagingstep4);
			 m_data_collectstress5[jj] = m_data_collectstress5[jj] + m_data_stress5[jj]/double(m_averagingstep4);
			m_data_collectvirial5[jj] = m_data_collectvirial5[jj] + m_data_virial5[jj]/double(m_averagingstep4);}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = m_data_collectD10[jj] + m_data_count10[jj]/double(m_averagingstep4);
			 m_data_collectT10[jj] = m_data_collectT10[jj] + m_data_temp10[jj]/double(m_averagingstep4);
			m_data_collectV10[jj] = m_data_collectV10[jj] + m_data_vel10[jj]/double(m_averagingstep4);
			m_data_collectstress10[jj] = m_data_collectstress10[jj] + m_data_stress10[jj]/double(m_averagingstep4);
			m_data_collectvirial10[jj] = m_data_collectvirial10[jj] + m_data_virial10[jj]/double(m_averagingstep4);}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = m_data_collectD50[jj] + m_data_count50[jj]/double(m_averagingstep4);
			 m_data_collectT50[jj] = m_data_collectT50[jj] + m_data_temp50[jj]/double(m_averagingstep4);
			 m_data_collectV50[jj] = m_data_collectV50[jj] + m_data_vel50[jj]/double(m_averagingstep4);}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = m_data_collectD100[jj] + m_data_count100[jj]/double(m_averagingstep4);
			 m_data_collectT100[jj] = m_data_collectT100[jj] + m_data_temp100[jj]/double(m_averagingstep4);
			m_data_collectV100[jj] = m_data_collectV100[jj] + m_data_vel100[jj]/double(m_averagingstep4);
			m_data_collectstress100[jj] = m_data_collectstress100[jj] + m_data_stress100[jj]/double(m_averagingstep4);
			m_data_collectvirial100[jj] = m_data_collectvirial100[jj] + m_data_virial100[jj]/double(m_averagingstep4);}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = m_data_collectD500[jj] + m_data_count500[jj]/double(m_averagingstep4);
			 m_data_collectT500[jj] = m_data_collectT500[jj] + m_data_temp500[jj]/double(m_averagingstep4);
			m_data_collectV500[jj] = m_data_collectV500[jj] + m_data_vel500[jj]/double(m_averagingstep4);}
		}

		if (d_ntimestep > m_startstep5 && d_ntimestep < m_endstep5+1)
		{ 
			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = m_data_collectD5[jj] + m_data_count5[jj]/double(m_averagingstep5);
			 m_data_collectT5[jj] = m_data_collectT5[jj] + m_data_temp5[jj]/double(m_averagingstep5);
			 m_data_collectV5[jj] = m_data_collectV5[jj] + m_data_vel5[jj]/double(m_averagingstep5);
			 m_data_collectstress5[jj] = m_data_collectstress5[jj] + m_data_stress5[jj]/double(m_averagingstep5);
			m_data_collectvirial5[jj] = m_data_collectvirial5[jj] + m_data_virial5[jj]/double(m_averagingstep5);}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = m_data_collectD10[jj] + m_data_count10[jj]/double(m_averagingstep5);
			 m_data_collectT10[jj] = m_data_collectT10[jj] + m_data_temp10[jj]/double(m_averagingstep5);
			m_data_collectV10[jj] = m_data_collectV10[jj] + m_data_vel10[jj]/double(m_averagingstep5);
			m_data_collectstress10[jj] = m_data_collectstress10[jj] + m_data_stress10[jj]/double(m_averagingstep5);
			m_data_collectvirial10[jj] = m_data_collectvirial10[jj] + m_data_virial10[jj]/double(m_averagingstep5);}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = m_data_collectD50[jj] + m_data_count50[jj]/double(m_averagingstep5);
			 m_data_collectT50[jj] = m_data_collectT50[jj] + m_data_temp50[jj]/double(m_averagingstep5);
			 m_data_collectV50[jj] = m_data_collectV50[jj] + m_data_vel50[jj]/double(m_averagingstep5);}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = m_data_collectD100[jj] + m_data_count100[jj]/double(m_averagingstep5);
			 m_data_collectT100[jj] = m_data_collectT100[jj] + m_data_temp100[jj]/double(m_averagingstep5);
			m_data_collectV100[jj] = m_data_collectV100[jj] + m_data_vel100[jj]/double(m_averagingstep5);
			m_data_collectstress100[jj] = m_data_collectstress100[jj] + m_data_stress100[jj]/double(m_averagingstep5);
			m_data_collectvirial100[jj] = m_data_collectvirial100[jj] + m_data_virial100[jj]/double(m_averagingstep5);}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = m_data_collectD500[jj] + m_data_count500[jj]/double(m_averagingstep5);
			 m_data_collectT500[jj] = m_data_collectT500[jj] + m_data_temp500[jj]/double(m_averagingstep5);
			m_data_collectV500[jj] = m_data_collectV500[jj] + m_data_vel500[jj]/double(m_averagingstep5);}
		}


		CString temp;
		
		// File Output
		if (d_ntimestep == double(m_endstep1+1.0))
		{
			m_ygridfile.Open(folder+"CollectD_5bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_5bin_1st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD5[ii], 18, collectD_5bin_1st);
					//m_ygridfile.WriteString(collectD_5bin_1st);
					temp.Format("%lf",m_data_collectD5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			
			m_ygridfile.Open(folder+"CollectT_5bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_5bin_1st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT5[ii], 18, collectT_5bin_1st);
					//m_ygridfile.WriteString(collectT_5bin_1st);
					temp.Format("%lf",m_data_collectT5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_5bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_5bin_1st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV5[ii], 18, collectV_5bin_1st);
					//m_ygridfile.WriteString(collectV_5bin_1st);
					temp.Format("%.15e",m_data_collectV5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_5bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_5bin_1st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress5[ii], 18, collectstress_5bin_1st);
					//m_ygridfile.WriteString(collectstress_5bin_1st);
					temp.Format("%.15e",m_data_collectstress5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			
			m_ygridfile.Open(folder+"Collectvirial_5bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_5bin_1st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial5[ii], 18, collectvirial_5bin_1st);
					//m_ygridfile.WriteString(collectvirial_5bin_1st);
					temp.Format("%.15e",m_data_collectvirial5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_10bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_10bin_1st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress10[ii], 18, collectstress_10bin_1st);
					//m_ygridfile.WriteString(collectstress_10bin_1st);
					temp.Format("%.15e",m_data_collectstress10[ii]);	m_ygridfile.WriteString(temp);				
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_10bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_10bin_1st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial10[ii], 18, collectvirial_10bin_1st);
					//m_ygridfile.WriteString(collectvirial_10bin_1st);
					temp.Format("%.15e",m_data_collectvirial10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_100bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_100bin_1st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress100[ii], 18, collectstress_100bin_1st);
					//m_ygridfile.WriteString(collectstress_100bin_1st);
					temp.Format("%.15e",m_data_collectstress100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_100bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_100bin_1st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial100[ii], 18, collectvirial_100bin_1st);
					//m_ygridfile.WriteString(collectvirial_100bin_1st);
					temp.Format("%.15e",m_data_collectvirial100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();





			m_ygridfile.Open(folder+"CollectD_10bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_10bin_1st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD10[ii], 18, collectD_10bin_1st);
					//m_ygridfile.WriteString(collectD_10bin_1st);
					temp.Format("%lf",m_data_collectD10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_10bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_10bin_1st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT10[ii], 18, collectT_10bin_1st);
					//m_ygridfile.WriteString(collectT_10bin_1st);
					temp.Format("%lf",m_data_collectT10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_10bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_10bin_1st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV10[ii], 18, collectV_10bin_1st);
					//m_ygridfile.WriteString(collectV_10bin_1st);
					temp.Format("%.15e",m_data_collectV10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_50bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_50bin_1st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD50[ii], 18, collectD_50bin_1st);
					//m_ygridfile.WriteString(collectD_50bin_1st);
					temp.Format("%lf",m_data_collectD50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_50bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_50bin_1st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT50[ii], 18, collectT_50bin_1st);
					//m_ygridfile.WriteString(collectT_50bin_1st);
					temp.Format("%lf",m_data_collectT50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_50bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_50bin_1st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV50[ii], 18, collectV_50bin_1st);
					//m_ygridfile.WriteString(collectV_50bin_1st);
					temp.Format("%.15e",m_data_collectV50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_100bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_100bin_1st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD100[ii], 18, collectD_100bin_1st);
					//m_ygridfile.WriteString(collectD_100bin_1st);
					temp.Format("%lf",m_data_collectD100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_100bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_100bin_1st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT100[ii], 18, collectT_100bin_1st);
					//m_ygridfile.WriteString(collectT_100bin_1st);
					temp.Format("%lf",m_data_collectT100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_100bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_100bin_1st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV100[ii], 18, collectV_100bin_1st);
					//m_ygridfile.WriteString(collectV_100bin_1st);
					temp.Format("%.15e",m_data_collectV100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_500bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_500bin_1st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD500[ii], 18, collectD_500bin_1st);
					//m_ygridfile.WriteString(collectD_500bin_1st);
					temp.Format("%lf",m_data_collectD500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_500bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_500bin_1st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT500[ii], 18, collectT_500bin_1st);
					//m_ygridfile.WriteString(collectT_500bin_1st);
					temp.Format("%lf",m_data_collectT500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_500bin_1st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_500bin_1st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV500[ii], 18, collectV_500bin_1st);
					//m_ygridfile.WriteString(collectV_500bin_1st);
					temp.Format("%.15e",m_data_collectV500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			

			for (jj = 0; jj < 5; jj++)
			{m_data_collectstress5[jj] = 0.0; m_data_collectvirial5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectstress10[jj] = 0.0; m_data_collectvirial10[jj] = 0.0;}

			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = 0.0; m_data_collectT5[jj] = 0.0; m_data_collectV5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = 0.0; m_data_collectT10[jj] = 0.0; m_data_collectV10[jj] = 0.0;}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = 0.0; m_data_collectT50[jj] = 0.0; m_data_collectV50[jj] = 0.0;}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = 0.0; m_data_collectT100[jj] = 0.0; m_data_collectV100[jj] = 0.0;
			m_data_collectstress100[jj] = 0.0; m_data_collectvirial100[jj] = 0.0;}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = 0.0; m_data_collectT500[jj] = 0.0; m_data_collectV500[jj] = 0.0;}
		}



		// File Output2
		if (d_ntimestep == double(m_endstep2+1.0))
		{
			m_ygridfile.Open(folder+"CollectD_5bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_5bin_2st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD5[ii], 18, collectD_5bin_2st);
					//m_ygridfile.WriteString(collectD_5bin_2st);
					temp.Format("%lf",m_data_collectD5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_5bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_5bin_2st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT5[ii], 18, collectT_5bin_2st);
					//m_ygridfile.WriteString(collectT_5bin_2st);
					temp.Format("%lf",m_data_collectT5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_5bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_5bin_2st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV5[ii], 18, collectV_5bin_2st);
					//m_ygridfile.WriteString(collectV_5bin_2st);
					temp.Format("%.15e",m_data_collectV5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_10bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_10bin_2st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress10[ii], 18, collectstress_10bin_2st);
					//m_ygridfile.WriteString(collectstress_10bin_2st);
					temp.Format("%.15e",m_data_collectstress10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_10bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_10bin_2st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial10[ii], 18, collectvirial_10bin_2st);
					//m_ygridfile.WriteString(collectvirial_10bin_2st);
					temp.Format("%.15e",m_data_collectvirial10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_100bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_100bin_2st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress100[ii], 18, collectstress_100bin_2st);
					//m_ygridfile.WriteString(collectstress_100bin_2st);
					temp.Format("%.15e",m_data_collectstress100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_100bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_100bin_2st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial100[ii], 18, collectvirial_100bin_2st);
					//m_ygridfile.WriteString(collectvirial_100bin_2st);
					temp.Format("%.15e",m_data_collectvirial100[ii]);	m_ygridfile.WriteString(temp);	
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_5bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_5bin_2st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress5[ii], 18, collectstress_5bin_2st);
					//m_ygridfile.WriteString(collectstress_5bin_2st);
					temp.Format("%.15e",m_data_collectstress5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_5bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_5bin_2st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial5[ii], 18, collectvirial_5bin_2st);
					//m_ygridfile.WriteString(collectvirial_5bin_2st);
					temp.Format("%.15e",m_data_collectvirial5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_10bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_10bin_2st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD10[ii], 18, collectD_10bin_2st);
					//m_ygridfile.WriteString(collectD_10bin_2st);
					temp.Format("%lf",m_data_collectD10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_10bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_10bin_2st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT10[ii], 18, collectT_10bin_2st);
					//m_ygridfile.WriteString(collectT_10bin_2st);
					temp.Format("%lf",m_data_collectT10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_10bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_10bin_2st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV10[ii], 18, collectV_10bin_2st);
					//m_ygridfile.WriteString(collectV_10bin_2st);
					temp.Format("%.15e",m_data_collectV10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_50bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_50bin_2st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD50[ii], 18, collectD_50bin_2st);
					//m_ygridfile.WriteString(collectD_50bin_2st);
					temp.Format("%lf",m_data_collectD50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_50bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_50bin_2st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT50[ii], 18, collectT_50bin_2st);
					//m_ygridfile.WriteString(collectT_50bin_2st);
					temp.Format("%lf",m_data_collectT50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_50bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_50bin_2st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV50[ii], 18, collectV_50bin_2st);
					//m_ygridfile.WriteString(collectV_50bin_2st);
					temp.Format("%.15e",m_data_collectV50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_100bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_100bin_2st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD100[ii], 18, collectD_100bin_2st);
					//m_ygridfile.WriteString(collectD_100bin_2st);
					temp.Format("%lf",m_data_collectD100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_100bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_100bin_2st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT100[ii], 18, collectT_100bin_2st);
					//m_ygridfile.WriteString(collectT_100bin_2st);
					temp.Format("%lf",m_data_collectT100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_100bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_100bin_2st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV100[ii], 18, collectV_100bin_2st);
					//m_ygridfile.WriteString(collectV_100bin_2st);
					temp.Format("%.15e",m_data_collectV100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_500bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_500bin_2st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD500[ii], 18, collectD_500bin_2st);
					//m_ygridfile.WriteString(collectD_500bin_2st);
					temp.Format("%lf",m_data_collectD500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_500bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_500bin_2st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT500[ii], 18, collectT_500bin_2st);
					//m_ygridfile.WriteString(collectT_500bin_2st);
					temp.Format("%lf",m_data_collectT500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_500bin_2st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_500bin_2st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV500[ii], 18, collectV_500bin_2st);
					//m_ygridfile.WriteString(collectV_500bin_2st);
					temp.Format("%.15e",m_data_collectV500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			for (jj = 0; jj < 5; jj++)
			{m_data_collectstress5[jj] = 0.0; m_data_collectvirial5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectstress10[jj] = 0.0; m_data_collectvirial10[jj] = 0.0;}

			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = 0.0; m_data_collectT5[jj] = 0.0; m_data_collectV5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = 0.0; m_data_collectT10[jj] = 0.0; m_data_collectV10[jj] = 0.0;}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = 0.0; m_data_collectT50[jj] = 0.0; m_data_collectV50[jj] = 0.0;}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = 0.0; m_data_collectT100[jj] = 0.0; m_data_collectV100[jj] = 0.0;
			m_data_collectstress100[jj] = 0.0; m_data_collectvirial100[jj] = 0.0;}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = 0.0; m_data_collectT500[jj] = 0.0; m_data_collectV500[jj] = 0.0;}
		}


		// File Output 3
		if (d_ntimestep == double(m_endstep3+1.0))
		{
			m_ygridfile.Open(folder+"CollectD_5bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_5bin_3rd[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD5[ii], 18, collectD_5bin_3rd);
					//m_ygridfile.WriteString(collectD_5bin_3rd);
					temp.Format("%lf",m_data_collectD5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_5bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_5bin_3rd[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT5[ii], 18, collectT_5bin_3rd);
					//m_ygridfile.WriteString(collectT_5bin_3rd);
					temp.Format("%lf",m_data_collectT5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_5bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_5bin_3rd[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV5[ii], 18, collectV_5bin_3rd);
					//m_ygridfile.WriteString(collectV_5bin_3rd);
					temp.Format("%.15e",m_data_collectV5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_5bin_3st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_5bin_3st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress5[ii], 18, collectstress_5bin_3st);
					//m_ygridfile.WriteString(collectstress_5bin_3st);
					temp.Format("%.15e",m_data_collectstress5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_5bin_3st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_5bin_3st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial5[ii], 18, collectvirial_5bin_3st);
					//m_ygridfile.WriteString(collectvirial_5bin_3st);
					temp.Format("%.15e",m_data_collectvirial5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_10bin_3st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_10bin_3st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress10[ii], 18, collectstress_10bin_3st);
					//m_ygridfile.WriteString(collectstress_10bin_3st);
					temp.Format("%.15e",m_data_collectstress10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_10bin_3st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_10bin_3st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial10[ii], 18, collectvirial_10bin_3st);
					//m_ygridfile.WriteString(collectvirial_10bin_3st);
					temp.Format("%.15e",m_data_collectvirial10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_100bin_3st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_100bin_3st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress100[ii], 18, collectstress_100bin_3st);
					//m_ygridfile.WriteString(collectstress_100bin_3st);
					temp.Format("%.15e",m_data_collectstress100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_100bin_3st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_100bin_3st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial100[ii], 18, collectvirial_100bin_3st);
					//m_ygridfile.WriteString(collectvirial_100bin_3st);
					temp.Format("%.15e",m_data_collectvirial100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_10bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_10bin_3rd[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD10[ii], 18, collectD_10bin_3rd);
					//m_ygridfile.WriteString(collectD_10bin_3rd);
					temp.Format("%lf",m_data_collectD10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_10bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_10bin_3rd[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT10[ii], 18, collectT_10bin_3rd);
					//m_ygridfile.WriteString(collectT_10bin_3rd);
					temp.Format("%lf",m_data_collectT10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_10bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_10bin_3rd[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV10[ii], 18, collectV_10bin_3rd);
					//m_ygridfile.WriteString(collectV_10bin_3rd);
					temp.Format("%.15e",m_data_collectV10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_50bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_50bin_3rd[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD50[ii], 18, collectD_50bin_3rd);
					//m_ygridfile.WriteString(collectD_50bin_3rd);
					temp.Format("%lf",m_data_collectD50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_50bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_50bin_3rd[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT50[ii], 18, collectT_50bin_3rd);
					//m_ygridfile.WriteString(collectT_50bin_3rd);
					temp.Format("%lf",m_data_collectT50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectV_50bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_50bin_3rd[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV50[ii], 18, collectV_50bin_3rd);
					//m_ygridfile.WriteString(collectV_50bin_3rd);
					temp.Format("%.15e",m_data_collectV50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_100bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_100bin_3rd[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD100[ii], 18, collectD_100bin_3rd);
					//m_ygridfile.WriteString(collectD_100bin_3rd);
					temp.Format("%lf",m_data_collectD100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_100bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_100bin_3rd[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT100[ii], 18, collectT_100bin_3rd);
					//m_ygridfile.WriteString(collectT_100bin_3rd);
					temp.Format("%lf",m_data_collectT100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_100bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_100bin_3rd[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV100[ii], 18, collectV_100bin_3rd);
					//m_ygridfile.WriteString(collectV_100bin_3rd);
					temp.Format("%.15e",m_data_collectV100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_500bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_500bin_3rd[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD500[ii], 18, collectD_500bin_3rd);
					//m_ygridfile.WriteString(collectD_500bin_3rd);
					temp.Format("%lf",m_data_collectD500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_500bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_500bin_3rd[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT500[ii], 18, collectT_500bin_3rd);
					//m_ygridfile.WriteString(collectT_500bin_3rd);
					temp.Format("%lf",m_data_collectT500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_500bin_3rd.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_500bin_3rd[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV500[ii], 18, collectV_500bin_3rd);
					//m_ygridfile.WriteString(collectV_500bin_3rd);
					temp.Format("%.15e",m_data_collectV500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			for (jj = 0; jj < 5; jj++)
			{m_data_collectstress5[jj] = 0.0; m_data_collectvirial5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectstress10[jj] = 0.0; m_data_collectvirial10[jj] = 0.0;}

			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = 0.0; m_data_collectT5[jj] = 0.0; m_data_collectV5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = 0.0; m_data_collectT10[jj] = 0.0; m_data_collectV10[jj] = 0.0;}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = 0.0; m_data_collectT50[jj] = 0.0; m_data_collectV50[jj] = 0.0;}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = 0.0; m_data_collectT100[jj] = 0.0; m_data_collectV100[jj] = 0.0;
			m_data_collectstress100[jj] = 0.0; m_data_collectvirial100[jj] = 0.0;}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = 0.0; m_data_collectT500[jj] = 0.0; m_data_collectV500[jj] = 0.0;}
		}


		// File Output 4
		if (d_ntimestep == double(m_endstep4+1.0))
		{
			m_ygridfile.Open(folder+"CollectD_5bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_5bin_4st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD5[ii], 18, collectD_5bin_4st);
					//m_ygridfile.WriteString(collectD_5bin_4st);
					temp.Format("%lf",m_data_collectD5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_5bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_5bin_4st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT5[ii], 18, collectT_5bin_4st);
					//m_ygridfile.WriteString(collectT_5bin_4st);
					temp.Format("%lf",m_data_collectT5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_5bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_5bin_4st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV5[ii], 18, collectV_5bin_4st);
					//m_ygridfile.WriteString(collectV_5bin_4st);
					temp.Format("%.15e",m_data_collectV5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_5bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_5bin_4st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress5[ii], 18, collectstress_5bin_4st);
					//m_ygridfile.WriteString(collectstress_5bin_4st);
					temp.Format("%.15e",m_data_collectstress5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_5bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_5bin_4st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial5[ii], 18, collectvirial_5bin_4st);
					//m_ygridfile.WriteString(collectvirial_5bin_4st);
					temp.Format("%.15e",m_data_collectvirial5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_10bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_10bin_4st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress10[ii], 18, collectstress_10bin_4st);
					//m_ygridfile.WriteString(collectstress_10bin_4st);
					temp.Format("%.15e",m_data_collectstress10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_10bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_10bin_4st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial10[ii], 18, collectvirial_10bin_4st);
					//m_ygridfile.WriteString(collectvirial_10bin_4st);
					temp.Format("%.15e",m_data_collectvirial10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_100bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_100bin_4st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress100[ii], 18, collectstress_100bin_4st);
					//m_ygridfile.WriteString(collectstress_100bin_4st);
					temp.Format("%.15e",m_data_collectstress100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_100bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_100bin_4st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial100[ii], 18, collectvirial_100bin_4st);
					//m_ygridfile.WriteString(collectvirial_100bin_4st);
					temp.Format("%.15e",m_data_collectvirial100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();




			m_ygridfile.Open(folder+"CollectD_10bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_10bin_4st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD10[ii], 18, collectD_10bin_4st);
					//m_ygridfile.WriteString(collectD_10bin_4st);
					temp.Format("%lf",m_data_collectD10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_10bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_10bin_4st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT10[ii], 18, collectT_10bin_4st);
					//m_ygridfile.WriteString(collectT_10bin_4st);
					temp.Format("%lf",m_data_collectT10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_10bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_10bin_4st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV10[ii], 18, collectV_10bin_4st);
					//m_ygridfile.WriteString(collectV_10bin_4st);
					temp.Format("%.15e",m_data_collectV10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_50bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_50bin_4st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD50[ii], 18, collectD_50bin_4st);
					//m_ygridfile.WriteString(collectD_50bin_4st);
					temp.Format("%lf",m_data_collectD50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_50bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_50bin_4st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT50[ii], 18, collectT_50bin_4st);
					//m_ygridfile.WriteString(collectT_50bin_4st);
					temp.Format("%lf",m_data_collectT50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_50bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_50bin_4st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV50[ii], 18, collectV_50bin_4st);
					//m_ygridfile.WriteString(collectV_50bin_4st);
					temp.Format("%.15e",m_data_collectV50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_100bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_100bin_4st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD100[ii], 18, collectD_100bin_4st);
					//m_ygridfile.WriteString(collectD_100bin_4st);
					temp.Format("%lf",m_data_collectD100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_100bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_100bin_4st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT100[ii], 18, collectT_100bin_4st);
					//m_ygridfile.WriteString(collectT_100bin_4st);
					temp.Format("%lf",m_data_collectT100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_100bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_100bin_4st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV100[ii], 18, collectV_100bin_4st);
					//m_ygridfile.WriteString(collectV_100bin_4st);
					temp.Format("%.15e",m_data_collectV100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_500bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_500bin_4st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD500[ii], 18, collectD_500bin_4st);
					//m_ygridfile.WriteString(collectD_500bin_4st);
					temp.Format("%lf",m_data_collectD500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_500bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_500bin_4st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT500[ii], 18, collectT_500bin_4st);
					//m_ygridfile.WriteString(collectT_500bin_4st);
					temp.Format("%lf",m_data_collectT500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_500bin_4st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_500bin_4st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV500[ii], 18, collectV_500bin_4st);
					//m_ygridfile.WriteString(collectV_500bin_4st);
					temp.Format("%.15e",m_data_collectV500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			for (jj = 0; jj < 5; jj++)
			{m_data_collectstress5[jj] = 0.0; m_data_collectvirial5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectstress10[jj] = 0.0; m_data_collectvirial10[jj] = 0.0;}


			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = 0.0; m_data_collectT5[jj] = 0.0; m_data_collectV5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = 0.0; m_data_collectT10[jj] = 0.0; m_data_collectV10[jj] = 0.0;}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = 0.0; m_data_collectT50[jj] = 0.0; m_data_collectV50[jj] = 0.0;}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = 0.0; m_data_collectT100[jj] = 0.0; m_data_collectV100[jj] = 0.0;
			m_data_collectstress100[jj] = 0.0; m_data_collectvirial100[jj] = 0.0;}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = 0.0; m_data_collectT500[jj] = 0.0; m_data_collectV500[jj] = 0.0;}
		}



		// File Output 5
		if (d_ntimestep == double(m_endstep5+1.0))
		{
			m_ygridfile.Open(folder+"CollectD_5bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_5bin_5st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD5[ii], 18, collectD_5bin_5st);
					//m_ygridfile.WriteString(collectD_5bin_5st);
					temp.Format("%lf",m_data_collectD5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_5bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_5bin_5st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT5[ii], 18, collectT_5bin_5st);
					//m_ygridfile.WriteString(collectT_5bin_5st);
					temp.Format("%lf",m_data_collectT5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_5bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_5bin_5st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV5[ii], 18, collectV_5bin_5st);
					//m_ygridfile.WriteString(collectV_5bin_5st);
					temp.Format("%.15e",m_data_collectV5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_5bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_5bin_5st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress5[ii], 18, collectstress_5bin_5st);
					//m_ygridfile.WriteString(collectstress_5bin_5st);
					temp.Format("%.15e",m_data_collectstress5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_5bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_5bin_5st[20];
			for(ii = 0; ii<5 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial5[ii], 18, collectvirial_5bin_5st);
					//m_ygridfile.WriteString(collectvirial_5bin_5st);
					temp.Format("%.15e",m_data_collectvirial5[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectstress_10bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_10bin_5st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress10[ii], 18, collectstress_10bin_5st);
					//m_ygridfile.WriteString(collectstress_10bin_5st);
					temp.Format("%.15e",m_data_collectstress10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_10bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_10bin_5st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial10[ii], 18, collectvirial_10bin_5st);
					//m_ygridfile.WriteString(collectvirial_10bin_5st);
					temp.Format("%.15e",m_data_collectvirial10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"Collectstress_100bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectstress_100bin_5st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectstress100[ii], 18, collectstress_100bin_5st);
					//m_ygridfile.WriteString(collectstress_100bin_5st);
					temp.Format("%.15e",m_data_collectstress100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"Collectvirial_100bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectvirial_100bin_5st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectvirial100[ii], 18, collectvirial_100bin_5st);
					//m_ygridfile.WriteString(collectvirial_100bin_5st);
					temp.Format("%.15e",m_data_collectvirial100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();



			m_ygridfile.Open(folder+"CollectD_10bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_10bin_5st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD10[ii], 18, collectD_10bin_5st);
					//m_ygridfile.WriteString(collectD_10bin_5st);
					temp.Format("%lf",m_data_collectD10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_10bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_10bin_5st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT10[ii], 18, collectT_10bin_5st);
					//m_ygridfile.WriteString(collectT_10bin_5st);
					temp.Format("%lf",m_data_collectT10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_10bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_10bin_5st[20];
			for(ii = 0; ii<10 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV10[ii], 18, collectV_10bin_5st);
					//m_ygridfile.WriteString(collectV_10bin_5st);
					temp.Format("%.15e",m_data_collectV10[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_50bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_50bin_5st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD50[ii], 18, collectD_50bin_5st);
					//m_ygridfile.WriteString(collectD_50bin_5st);
					temp.Format("%lf",m_data_collectD50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_50bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_50bin_5st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT50[ii], 18, collectT_50bin_5st);
					//m_ygridfile.WriteString(collectT_50bin_5st);
					temp.Format("%lf",m_data_collectT50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_50bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_50bin_5st[20];
			for(ii = 0; ii<40 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV50[ii], 18, collectV_50bin_5st);
					//m_ygridfile.WriteString(collectV_50bin_5st);
					temp.Format("%.15e",m_data_collectV50[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			m_ygridfile.Open(folder+"CollectD_100bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_100bin_5st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD100[ii], 18, collectD_100bin_5st);
					//m_ygridfile.WriteString(collectD_100bin_5st);
					temp.Format("%lf",m_data_collectD100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_100bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_100bin_5st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT100[ii], 18, collectT_100bin_5st);
					//m_ygridfile.WriteString(collectT_100bin_5st);
					temp.Format("%lf",m_data_collectT100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_100bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_100bin_5st[20];
			for(ii = 0; ii<100 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV100[ii], 18, collectV_100bin_5st);
					//m_ygridfile.WriteString(collectV_100bin_5st);
					temp.Format("%.15e",m_data_collectV100[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();




			m_ygridfile.Open(folder+"CollectD_500bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectD_500bin_5st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectD500[ii], 18, collectD_500bin_5st);
					//m_ygridfile.WriteString(collectD_500bin_5st);
					temp.Format("%lf",m_data_collectD500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();
			m_ygridfile.Open(folder+"CollectT_500bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectT_500bin_5st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectT500[ii], 18, collectT_500bin_5st);
					//m_ygridfile.WriteString(collectT_500bin_5st);
					temp.Format("%lf",m_data_collectT500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();

			m_ygridfile.Open(folder+"CollectV_500bin_5st.txt",CFile::modeCreate | CFile::modeWrite);
			//char collectV_500bin_5st[20];
			for(ii = 0; ii<25 ; ii++) //Write x,y,vx,vy
			{
					//_gcvt(m_data_collectV500[ii], 18, collectV_500bin_5st);
					//m_ygridfile.WriteString(collectV_500bin_5st);
					temp.Format("%.15e",m_data_collectV500[ii]);	m_ygridfile.WriteString(temp);
					m_ygridfile.WriteString("\n"); 
			}
			m_ygridfile.Close();


			for (jj = 0; jj < 5; jj++)
			{m_data_collectstress5[jj] = 0.0; m_data_collectvirial5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectstress10[jj] = 0.0; m_data_collectvirial10[jj] = 0.0;}

			for (jj = 0; jj < 5; jj++)
			{m_data_collectD5[jj] = 0.0; m_data_collectT5[jj] = 0.0; m_data_collectV5[jj] = 0.0;}
			for (jj = 0; jj < 10; jj++)
			{m_data_collectD10[jj] = 0.0; m_data_collectT10[jj] = 0.0; m_data_collectV10[jj] = 0.0;}
			for (jj = 0; jj < 40; jj++)
			{m_data_collectD50[jj] = 0.0; m_data_collectT50[jj] = 0.0; m_data_collectV50[jj] = 0.0;}
			for (jj = 0; jj < 100; jj++)
			{m_data_collectD100[jj] = 0.0; m_data_collectT100[jj] = 0.0; m_data_collectV100[jj] = 0.0;
			m_data_collectstress100[jj] = 0.0; m_data_collectvirial100[jj] = 0.0;}
			for (jj = 0; jj < 25; jj++)
			{m_data_collectD500[jj] = 0.0; m_data_collectT500[jj] = 0.0; m_data_collectV500[jj] = 0.0;}

									d_ntimestep = 0;

									m_heatsumx = 0.0;
									m_heatsumy = 0.0;
									m_heatsumz = 0.0;
									d_totalLoop += 1;
									TmpFileSave();

									_itoa(Nfolder, bufferfolder, 10 );
									folder = CString(bufferfolder);
									Nfolder++; //number of folder used for save file
		}


//Collect Data For long time averaging collection


}



void CMFCGLView::SetInitialCondition()
{


			m_nMs = 4;


			m_aM[0].m_fX =2.3;
			m_aM[0].m_fY =2.3;
			m_aM[0].m_fZ =2.3;


			m_aM[1].m_fX =2.65;
			m_aM[1].m_fY =2.65;
			m_aM[1].m_fZ =2.65;

			m_aM[2].m_fX =1.3;
			m_aM[2].m_fY =1.3;
			m_aM[2].m_fZ =1.3;


			m_aM[3].m_fX =1.65;
			m_aM[3].m_fY =1.65;
			m_aM[3].m_fZ =1.65;


//				m_nMs = 2;
//
//			m_aM[0].m_fX =2.0;
//			m_aM[0].m_fY =2.0;
//			m_aM[0].m_fZ =1.3;
//
//			m_aM[1].m_fX =2.5;
//			m_aM[1].m_fY =2.0;
//			m_aM[1].m_fZ =1.3;
//						

	/*

			m_nMs = 27;

			m_aM[0].m_fX =2.0;
			m_aM[0].m_fY =2.0;
			m_aM[0].m_fZ =1.3;

			m_aM[1].m_fX =2.5;
			m_aM[1].m_fY =2.0;
			m_aM[1].m_fZ =1.3;
						
			m_aM[2].m_fX =3.0;
			m_aM[2].m_fY =2.0;
			m_aM[2].m_fZ =1.3;

			m_aM[3].m_fX =2.0;
			m_aM[3].m_fY =2.5;
			m_aM[3].m_fZ =1.3;

			m_aM[4].m_fX =2.5;
			m_aM[4].m_fY =2.5;
			m_aM[4].m_fZ =1.3;
						
			m_aM[5].m_fX =3.0;
			m_aM[5].m_fY =2.5;
			m_aM[5].m_fZ =1.3;

			m_aM[6].m_fX =2.0;
			m_aM[6].m_fY =3.0;
			m_aM[6].m_fZ =1.3;

			m_aM[7].m_fX =2.5;
			m_aM[7].m_fY =3.0;
			m_aM[7].m_fZ =1.3;
						
			m_aM[8].m_fX =3.0;
			m_aM[8].m_fY =3.0;
			m_aM[8].m_fZ =1.3;

			m_aM[9].m_fX =2.0;
			m_aM[9].m_fY =2.0;
			m_aM[9].m_fZ =1.9;

			m_aM[10].m_fX =2.5;
			m_aM[10].m_fY =2.0;
			m_aM[10].m_fZ =1.9;
						
			m_aM[11].m_fX =3.0;
			m_aM[11].m_fY =2.0;
			m_aM[11].m_fZ =1.9;

			m_aM[12].m_fX =2.0;
			m_aM[12].m_fY =2.6;
			m_aM[12].m_fZ =1.9;

			m_aM[13].m_fX =2.5;
			m_aM[13].m_fY =2.6;
			m_aM[13].m_fZ =1.9;
						
			m_aM[14].m_fX =3.0;
			m_aM[14].m_fY =2.6;
			m_aM[14].m_fZ =1.9;


			m_aM[15].m_fX =2.0;
			m_aM[15].m_fY =3.2;
			m_aM[15].m_fZ =1.9;

			m_aM[16].m_fX =2.5;
			m_aM[16].m_fY =3.2;
			m_aM[16].m_fZ =1.9;
						
			m_aM[17].m_fX =3.0;
			m_aM[17].m_fY =3.2;
			m_aM[17].m_fZ =1.9;


			m_aM[18].m_fX =2.0;
			m_aM[18].m_fY =2.0;
			m_aM[18].m_fZ =2.3;

			m_aM[19].m_fX =2.5;
			m_aM[19].m_fY =2.0;
			m_aM[19].m_fZ =2.3;
						
			m_aM[20].m_fX =3.0;
			m_aM[20].m_fY =2.0;
			m_aM[20].m_fZ =2.3;

			m_aM[21].m_fX =2.0;
			m_aM[21].m_fY =2.6;
			m_aM[21].m_fZ =2.3;

			m_aM[22].m_fX =2.5;
			m_aM[22].m_fY =2.6;
			m_aM[22].m_fZ =2.3;
						
			m_aM[23].m_fX =3.0;
			m_aM[23].m_fY =2.6;
			m_aM[23].m_fZ =2.3;

			m_aM[24].m_fX =2.0;
			m_aM[24].m_fY =3.2;
			m_aM[24].m_fZ =2.3;

			m_aM[25].m_fX =2.5;
			m_aM[25].m_fY =3.2;
			m_aM[25].m_fZ =2.3;
						
			m_aM[26].m_fX =3.0;
			m_aM[26].m_fY =3.2;
			m_aM[26].m_fZ =2.3;
*/

			for (int ii = 0 ; ii < m_nMs; ii ++)
			{
				m_aM[ii].m_fX =m_aM[ii].m_fX-0.9;
				m_aM[ii].m_fY =m_aM[ii].m_fY-0.9;
				m_aM[ii].m_fZ =m_aM[ii].m_fZ-0.9;
			}
			




}








void CMFCGLView::OnMove()
{
		// Change mode to coordinate change
	if (m_nMode == 0)
	{	m_nMode = 1;	} 
	else
	{	m_nMode = 0;	}
	
}

void CMFCGLView::OnShowmolecule()
{

	// ShowMolecule // Not use now
	if (m_nShow_molecule==0)
	{		m_nShow_molecule=1;	}
	else
	{		m_nShow_molecule=0;	}
	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
	
}

void CMFCGLView::OnNewsim()
{
		// TODO: Add your command handler code here
	m_nMode = 1;
	m_nShow_molecule=0;
	m_xRotate=0.0;
	m_yRotate=0.0;
	m_zRotate=0.0;
	m_nWire =0;
	m_Gcount = 0;
	d_ntimestep = 0.0;
	
	m_nSp = 0;

		for(int ii = 0; ii<m_nMs ; ii++)
		{
			m_aM[ii].m_fVx = 0.0;
			m_aM[ii].m_fVy = 0.0;
			m_aM[ii].m_fVz = 0.0;

			m_aM[ii].m_fX =0.0;
			m_aM[ii].m_fY =0.0;
			m_aM[ii].m_fZ =0.0;

			m_aM[ii].m_fax = 0.0;
			m_aM[ii].m_fay = 0.0;
			m_aM[ii].m_faz = 0.0;


		}
		m_nMs = 0;

	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
}

void CMFCGLView::OnWire()
{
		// TODO: Add your command handler code here
	m_nWire = m_nWire + 1;
	if (m_nWire==4)
	{	m_nWire = 0;	}
	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
}

void CMFCGLView::OnRun()
{
		// Locate molecule on the cell border and find jump in acceleration - BHK
	int maxroop;
	if (m_time == 1) maxroop = m_nroop; // Number of loop per one Run()
	else maxroop = 1;
	int ii;

	//if (m_nMs == 0) return;
	//if (m_nwMs == 0) return;

	if (d_ntimestep > m_endstep5 + 1000.0)
	{
		d_ntimestep = 0.0;
		m_heatsumx = 0.0;
		m_heatsumy = 0.0;
		m_heatsumz = 0.0;
		d_totalLoop += 1;
		TmpFileSave();
	}


	start = clock();


	for (int roop = 0 ; roop <maxroop ; roop++)  //One time step calculation
	{
		d_ntimestep = d_ntimestep + 1.0; //for record n'th timestep
		//set to 0 in order to off file saving

		// First timestep at velocity Verlet


// 0.5 v setp for fluid 
		for(ii = 0; ii<m_nMs ; ii++)
		{


			//Get next positon (Prevent inf acc)
			if (m_aM[ii].m_fax > 2000 )	
			{m_aM[ii].m_fax = 0;m_aM[ii].m_fay = 0;m_aM[ii].m_fVx = 0.0;m_aM[ii].m_fVy = 0.0;m_aM[ii].m_fVz = 0.0;}//AfxMessageBox("error"); break;}
			if (m_aM[ii].m_fax < -2000 )	
			{m_aM[ii].m_fax = 0;m_aM[ii].m_fay = 0;m_aM[ii].m_fVx = 0.0;m_aM[ii].m_fVy = 0.0;m_aM[ii].m_fVz = 0.0;}//AfxMessageBox("error"); break;}

			m_aM[ii].m_fX =m_aM[ii].m_fX + m_fdt*m_aM[ii].m_fVx + 0.5*m_fdt*m_fdt*m_aM[ii].m_fax;  //Verlet (8b) M P. Allen
			m_aM[ii].m_fY =m_aM[ii].m_fY + m_fdt*m_aM[ii].m_fVy + 0.5*m_fdt*m_fdt*m_aM[ii].m_fay;
			m_aM[ii].m_fZ =m_aM[ii].m_fZ + m_fdt*m_aM[ii].m_fVz + 0.5*m_fdt*m_fdt*m_aM[ii].m_faz;

//			m_aM[ii].m_fX =m_aM[ii].m_fX + m_fdt*m_aM[ii].m_fVx;  //Verlet (8b) M P. Allen
//			m_aM[ii].m_fY =m_aM[ii].m_fY + m_fdt*m_aM[ii].m_fVy;
//			m_aM[ii].m_fZ =m_aM[ii].m_fZ + m_fdt*m_aM[ii].m_fVz;

			m_aM[ii].m_fVx = m_aM[ii].m_fVx + 0.5*m_fdt*m_aM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aM[ii].m_fVy = m_aM[ii].m_fVy + 0.5*m_fdt*m_aM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aM[ii].m_fVz = m_aM[ii].m_fVz + 0.5*m_fdt*m_aM[ii].m_faz;  //Verlet (8a) M P. Allen
			//End of t+dt timestep

			m_aM[ii].m_fax = 0.0;
			m_aM[ii].m_fay = 0.0;
			m_aM[ii].m_faz = 0.0;
			m_aM[ii].m_fPot = 0.0;
			m_aM[ii].m_fStress = 0.0;
			m_aM[ii].m_fVirialxy = 0.0;

			m_aM_OMP[ii].m_fax = 0.0;
			m_aM_OMP[ii].m_fay = 0.0;
			m_aM_OMP[ii].m_faz = 0.0;
			m_aM_OMP[ii].m_fPot = 0.0;
			m_aM_OMP[ii].m_fStress = 0.0;
			m_aM_OMP[ii].m_fVirialxy = 0.0;

		}

		//// 0.5 v setp for wall 
		for( ii = 0; ii<m_nwMs ; ii++)
		{

			m_aTwM[ii].m_fX =m_aTwM[ii].m_fX + m_fdt*m_aTwM[ii].m_fVx + 0.5*m_fdt*m_fdt*m_aTwM[ii].m_fax;  //Verlet (8b) M P. Allen
			m_aTwM[ii].m_fY =m_aTwM[ii].m_fY + m_fdt*m_aTwM[ii].m_fVy + 0.5*m_fdt*m_fdt*m_aTwM[ii].m_fay;
			m_aTwM[ii].m_fZ =m_aTwM[ii].m_fZ + m_fdt*m_aTwM[ii].m_fVz + 0.5*m_fdt*m_fdt*m_aTwM[ii].m_faz;

			m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx + 0.5*m_fdt*m_aTwM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy + 0.5*m_fdt*m_aTwM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz + 0.5*m_fdt*m_aTwM[ii].m_faz;  //Verlet (8a) M P. Allen
			//End of t+dt timestep

//			m_aTwM[ii].m_fX =m_aTwM[ii].m_fX + m_fdt*m_aTwM[ii].m_fVx;  //Verlet (8b) M P. Allen
//			m_aTwM[ii].m_fY =m_aTwM[ii].m_fY + m_fdt*m_aTwM[ii].m_fVy;
//			m_aTwM[ii].m_fZ =m_aTwM[ii].m_fZ + m_fdt*m_aTwM[ii].m_fVz;

			m_aTwM[ii].m_fax = 0.0;
			m_aTwM[ii].m_fay = 0.0;
			m_aTwM[ii].m_faz = 0.0;
			m_aTwM[ii].m_fPot = 0.0;
			m_aTwM[ii].m_fStress = 0.0;
			m_aTwM[ii].m_fVirialxy = 0.0;

			m_aTwM_OMP[ii].m_fax = 0.0;
			m_aTwM_OMP[ii].m_fay = 0.0;
			m_aTwM_OMP[ii].m_faz = 0.0;
			m_aTwM_OMP[ii].m_fPot = 0.0;
			m_aTwM_OMP[ii].m_fStress = 0.0;
			m_aTwM_OMP[ii].m_fVirialxy = 0.0;
			


			m_aBwM[ii].m_fX =m_aBwM[ii].m_fX + m_fdt*m_aBwM[ii].m_fVx + 0.5*m_fdt*m_fdt*m_aBwM[ii].m_fax;  //Verlet (8b) M P. Allen
			m_aBwM[ii].m_fY =m_aBwM[ii].m_fY + m_fdt*m_aBwM[ii].m_fVy + 0.5*m_fdt*m_fdt*m_aBwM[ii].m_fay;
			m_aBwM[ii].m_fZ =m_aBwM[ii].m_fZ + m_fdt*m_aBwM[ii].m_fVz + 0.5*m_fdt*m_fdt*m_aBwM[ii].m_faz;

			m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx + 0.5*m_fdt*m_aBwM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy + 0.5*m_fdt*m_aBwM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz + 0.5*m_fdt*m_aBwM[ii].m_faz;  //Verlet (8a) M P. Allen
			//End of t+dt timestep

//			m_aBwM[ii].m_fX =m_aBwM[ii].m_fX + m_fdt*m_aBwM[ii].m_fVx;  //Verlet (8b) M P. Allen
//			m_aBwM[ii].m_fY =m_aBwM[ii].m_fY + m_fdt*m_aBwM[ii].m_fVy;
//			m_aBwM[ii].m_fZ =m_aBwM[ii].m_fZ + m_fdt*m_aBwM[ii].m_fVz;


			m_aBwM[ii].m_fax = 0.0;
			m_aBwM[ii].m_fay = 0.0;
			m_aBwM[ii].m_faz = 0.0;
			m_aBwM[ii].m_fPot = 0.0;
			m_aBwM[ii].m_fStress = 0.0;
			m_aBwM[ii].m_fVirialxy = 0.0;
						
			m_aBwM_OMP[ii].m_fax = 0.0;
			m_aBwM_OMP[ii].m_fay = 0.0;
			m_aBwM_OMP[ii].m_faz = 0.0;
			m_aBwM_OMP[ii].m_fPot = 0.0;
			m_aBwM_OMP[ii].m_fStress = 0.0;
			m_aBwM_OMP[ii].m_fVirialxy = 0.0;

		}


		SortGridInit();
//Gas Molecule Interaction : Cell algorithm.
		GetCellInteraction();

									//Reserve for OPENMP	
									//int op;
									//for(op = 0; op<m_nMs ; op++) //this is for open mp, update outside OpenMP parallel loop
									//{
									//	m_aM[op].m_fax = m_aM_OMP[op].m_fax;
									//	m_aM[op].m_fay = m_aM_OMP[op].m_fay;
									//	m_aM[op].m_faz = m_aM_OMP[op].m_faz;
									//	m_aM[op].m_fPot = m_aM_OMP[op].m_fPot;
									//	m_aM[op].m_fStress = m_aM_OMP[op].m_fStress;
									//	m_aM[op].m_fVirialxy = m_aM_OMP[op].m_fVirialxy;
									//}


	    for( ii = 0; ii<m_nwMs ; ii++)
		{
		m_aTwM[ii].m_fax =m_aTwM[ii].m_fax + d_K_M*(-m_aTwM[ii].m_fX + m_aTwM[ii].m_fXo);
		m_aTwM[ii].m_fay =m_aTwM[ii].m_fay + d_K_M*(-m_aTwM[ii].m_fY + m_aTwM[ii].m_fYo);
		m_aTwM[ii].m_faz =m_aTwM[ii].m_faz + d_K_M*(-m_aTwM[ii].m_fZ + m_aTwM[ii].m_fZo);

		m_aBwM[ii].m_fax =m_aBwM[ii].m_fax + d_K_M*(-m_aBwM[ii].m_fX + m_aBwM[ii].m_fXo);
		m_aBwM[ii].m_fay =m_aBwM[ii].m_fay + d_K_M*(-m_aBwM[ii].m_fY + m_aBwM[ii].m_fYo);
		m_aBwM[ii].m_faz =m_aBwM[ii].m_faz + d_K_M*(-m_aBwM[ii].m_fZ + m_aBwM[ii].m_fZo);
		}

		//For pressure/shear driven flow;
		DriveFlow();  


//Gas Molecule Interaction : Cell algorithm. -END

	//! Important
	m_vec.fsignx = 0.0;	m_vec.fsigny = 0.0; m_vec.fsigny = 0.0;
	m_vec.accx = 0.0; m_vec.accy = 0.0; m_vec.accz = 0.0;
	m_vec.fForce = 0.0; m_vec.fPot = 0.0;


//Nose-Hoover
//	m_zeta = 0.0;
//	m_zeta_dot = 0.0;
//	m_tau = 0.0;
//	m_Q = 0.0;

	m_tau = 0.09622;  //0.1ps -> D.J Evans, and Brad Lee Holian paper 1985 

	if (m_T != 0)
	{
	m_Q = 3.0*double(m_nMs-1)*m_Ttarget*m_tau*m_tau;
	m_zeta_dot = (2.0*m_Kin*m_nMs-3.0*(m_nMs-1)*m_Ttarget)/m_Q;
	//m_zeta_dot = (2.0*m_Kin*m_nMs-2.0*(m_nMs-1)*m_Ttarget)/m_Q; //for exclude flow direction 1/2
	m_zeta = m_zeta + m_zeta_dot*m_fdt;
	}

//	m_T = 0.5*m_Kb*m_V/double(m_nMs);
//	m_Kin = 0.5*m_Kb*m_V/double(m_nMs);
//Nose-Hoover - END


	//Get velocity -> move to next position
	    for( ii = 0; ii<m_nMs ; ii++)
		{
			//Get velocity
			//V(n+3/2) = V(n+1/2) + ai Δt
			//V(n+1) = (V(n+3/2) +V(n+1/2))/2

			if (m_ndamping == 2)
			{
			//With Nose-Hoover
			m_aM[ii].m_fVx = m_aM[ii].m_fVx*(1.0-m_zeta*0.5*m_fdt) + 0.5*m_fdt*m_aM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aM[ii].m_fVy = m_aM[ii].m_fVy*(1.0-m_zeta*0.5*m_fdt) + 0.5*m_fdt*m_aM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aM[ii].m_fVz = m_aM[ii].m_fVz*(1.0-m_zeta*0.5*m_fdt) + 0.5*m_fdt*m_aM[ii].m_faz;  //Verlet (8a) M P. Allen
			}else
			{
			//Velocity at t+dt timestep
			m_aM[ii].m_fVx = m_aM[ii].m_fVx + 0.5*m_fdt*m_aM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aM[ii].m_fVy = m_aM[ii].m_fVy + 0.5*m_fdt*m_aM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aM[ii].m_fVz = m_aM[ii].m_fVz + 0.5*m_fdt*m_aM[ii].m_faz;  //Verlet (8a) M P. Allen
			}
			

		}

	    for( ii = 0; ii<m_nwMs ; ii++)
		{
			//V(n+3/2) = V(n+1/2) + ai Δt
			//V(n+1) = (V(n+3/2) +V(n+1/2))/2
			//Velocity at t+dt timestep
			m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx + 0.5*m_fdt*m_aTwM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy + 0.5*m_fdt*m_aTwM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz + 0.5*m_fdt*m_aTwM[ii].m_faz;  //Verlet (8a) M P. Allen

			m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx + 0.5*m_fdt*m_aBwM[ii].m_fax;  //m_fVxh = m_fVx(t+0.5dt) 
			m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy + 0.5*m_fdt*m_aBwM[ii].m_fay;  //Verlet (8a) M P. Allen
			m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz + 0.5*m_fdt*m_aBwM[ii].m_faz;  //Verlet (8a) M P. Allen
		}


		//Collecting velocity data
		CollectAndSaveFile(); //For file output.

		CollectData(); //HERE IS previous SET Not t+0.5dt V and t-Pot  FOR screen output

		//Uncomment here for heat flux calculation
		CollectEnergy(); //for heat flux calculation.
		
		
		//Bohung Here

		for( ii = 0; ii<m_nMs ; ii++)
		{

			//Damping -Reduce velocity
			if (m_ndamping == 1)
			{
				if (m_T != 0.0)	m_dampcoef = sqrt(m_Ttarget/m_T);

				m_aM[ii].m_fVx = m_aM[ii].m_fVx*m_dampcoef;
				m_aM[ii].m_fVy = m_aM[ii].m_fVy*m_dampcoef;
				m_aM[ii].m_fVz = m_aM[ii].m_fVz*m_dampcoef;
			}
			if (m_ndamping == 0)
			{
				m_aM[ii].m_fVx = m_aM[ii].m_fVx/m_dampcoef; 
				m_aM[ii].m_fVy = m_aM[ii].m_fVy/m_dampcoef;
				m_aM[ii].m_fVz = m_aM[ii].m_fVz/m_dampcoef;
			}
		}

  

		double walldampT,walldampB;
		walldampB = 1.0;
		walldampT = 1.0;
		//walldamp = 1.0;

		if (wallTempT != 0.0)  //never block
		{
			double tempy;
			tempy = m_fH+(m_fMdist/2.0);
			for( ii = 0; ii<m_nwMs ; ii++)
			{
				walldampT = 1.0;
				if (d_Twalltemp1 != 0.0) walldampT = sqrt(wallTempT/(d_Twalltemp1));
				if (m_aTwM[ii].m_fYo - tempy < 0.0) 
				{
					m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx*walldampT;
					m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy*walldampT;
					m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz*walldampT;
				}

				walldampT = 1.0;
				if (d_Twalltemp2 != 0.0) walldampT = sqrt(wallTempT/(d_Twalltemp2));
				if (fabs(m_aTwM[ii].m_fYo - tempy) < m_fMdist/3.0) 
				{
					m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx*walldampT;
					m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy*walldampT;
					m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz*walldampT;
				}

				walldampT = 1.0;					
				if (d_Twalltemp3 != 0.0) walldampT = sqrt(wallTempT/(d_Twalltemp3));
				if (m_aTwM[ii].m_fYo - tempy > 0.0)
				{
					m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx*walldampT;
					m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy*walldampT;
					m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz*walldampT;
				}
				
					
			}
		}else //temp == 0 이면 
		{
			for( ii = 0; ii<m_nwMs ; ii++)
			{
			m_aTwM[ii].m_fVx = 0.0;
			m_aTwM[ii].m_fVy = 0.0;
			m_aTwM[ii].m_fVz = 0.0;
			m_aTwM[ii].m_fax = 0.0;
			m_aTwM[ii].m_fay = 0.0;
			m_aTwM[ii].m_faz = 0.0;
			m_aTwM[ii].m_fX = m_aTwM[ii].m_fXo;
			m_aTwM[ii].m_fY = m_aTwM[ii].m_fYo;
			m_aTwM[ii].m_fZ = m_aTwM[ii].m_fZo;
			}
		}

		if (wallTempB != 0.0)  //never block
		{

			for( ii = 0; ii<m_nwMs ; ii++)
			{
				walldampB = 1.0;
				if (d_Bwalltemp1 != 0.0) walldampB = sqrt(wallTempB/(d_Bwalltemp1));
				if (m_aBwM[ii].m_fYo == 0.0) 
				{
					m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx*walldampB;
					m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy*walldampB;
					m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz*walldampB;
				}

				walldampB = 1.0;
				if (d_Bwalltemp2 != 0.0) walldampB = sqrt(wallTempB/(d_Bwalltemp2));
				if (m_aBwM[ii].m_fYo == -m_fMdist/2.0) 
				{
					m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx*walldampB;
					m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy*walldampB;
					m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz*walldampB;
				}

				walldampB = 1.0;					
				if (d_Bwalltemp3 != 0.0) walldampB = sqrt(wallTempB/(d_Bwalltemp3));
				if (m_aBwM[ii].m_fYo == -m_fMdist) 
				{
					m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx*walldampB;
					m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy*walldampB;
					m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz*walldampB;
				}
			}

		}else  //temp == 0 이면 
		{
			for( ii = 0; ii<m_nwMs ; ii++)
			{
			m_aBwM[ii].m_fVx = 0.0;
			m_aBwM[ii].m_fVy = 0.0;
			m_aBwM[ii].m_fVz = 0.0;
			m_aBwM[ii].m_fax = 0.0;
			m_aBwM[ii].m_fay = 0.0;
			m_aBwM[ii].m_faz = 0.0;
			m_aBwM[ii].m_fX = m_aBwM[ii].m_fXo;
			m_aBwM[ii].m_fY = m_aBwM[ii].m_fYo;
			m_aBwM[ii].m_fZ = m_aBwM[ii].m_fZo;
			}
		}



	
	}//for (int roop = 0 ; roop <maxroop ; roop++)  //One time step calculation -END



	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	

		m_Potsca1[m_Gcount] = 0.005*m_Kb*m_aM[0].m_fPot;
		m_Potsca2[m_Gcount] = 0.005*m_Kb*m_aM[1].m_fPot;
		m_Potsca3[m_Gcount] = 0.005*m_Kb*m_aM[2].m_fPot;
		m_Potsca4[m_Gcount] = 0.005*m_Kb*m_aM[3].m_fPot;

		m_Kinsca1[m_Gcount] = 0.005*m_Kb*1.0 * (pow(m_aM[0].m_fVx,2.0) + pow(m_aM[0].m_fVy,2.0));
		m_Kinsca2[m_Gcount] = 0.005*m_Kb*1.0 * (pow(m_aM[1].m_fVx,2.0) + pow(m_aM[1].m_fVy,2.0));
		m_Kinsca3[m_Gcount] = 0.005*m_Kb*1.0 * (pow(m_aM[2].m_fVx,2.0) + pow(m_aM[2].m_fVy,2.0));
		m_Kinsca4[m_Gcount] = 0.005*m_Kb*1.0 * (pow(m_aM[3].m_fVx,2.0) + pow(m_aM[3].m_fVy,2.0));

		m_Tsc   = m_T/100.0;
		m_Kinsc = m_Kin/200.0;
		m_Potsc = m_Pot/200.0;
		m_Totsc = m_Tot/200.0;

		m_Gcount = m_Gcount + 1;

		m_Tsca[m_Gcount] = m_Tsc;
		m_Kinsca[m_Gcount] = m_Kinsc;
		m_Potsca[m_Gcount] = m_Potsc;
		m_Totsca[m_Gcount] = m_Totsc;
		if (m_Gcount == 399) m_Gcount = 0;


	//DrawScene(); SwapBuffers(m_pDC->GetSafeHdc()); //Update the Screen

////MB////////
if (visual == 0)	DrawScene();
////MB///////

//	DrawScene();
	SwapBuffers(m_pDC->GetSafeHdc());
//	Invalidate(FALSE);
	//Change to swapbuffer()
	
}

void CMFCGLView::OnShear()
{
		// TODO: Add your command handler code here
	if (m_nSp != 2)
	{
	m_nSp = 2;
	m_gamma = 2.0*m_Umax/m_fH;
 //m_data_temp10[jj] = m_data_temp10[jj] + pow((m_aM[ii].m_fVx -m_gamma2*pBV - m_gamma*(m_aM[ii].m_fY-m_fH/2.0)),2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
	}
	else
	{ m_nSp = 0;
	  m_gamma = 0.0;
	}
}

void CMFCGLView::OnPeriodic()
{
		// TODO: Add your command handler code here
	if (m_nPeriodic == 0) 
	{	m_nPeriodic =1;	}
	else
	{	m_nPeriodic =0;	}
}

void CMFCGLView::OnPressure()
{
		// TODO: Add your command handler code here
	if (m_nSp != 1)
	{
		m_nSp = 1;
		m_gamma2 = 1;
	}
	else
	{
		m_nSp = 0;		
		m_gamma2 = 0;
	}
}

void CMFCGLView::OnInit()
{
		// TODO: Add your command handler code here
	//
	//CString temp;
	//double testnum = 0.00000000000012341235213412352135;
	//temp.Format("%.15e",testnum);
	//AfxMessageBox(temp);

	SetInitialCondition();
	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
}

void CMFCGLView::OnPropdlg()
{

		if (m_time == 1)
	{
		m_time = 0;
		KillTimer(1);
		m_dautocolor = 0.4;
	}
	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());

	//OnNewsim();

	//TODO: Add your command handler code here
	//m_Propdlg.m_nwMs = m_nwMs;
	m_Propdlg.m_nwMs = m_nwMs;
	m_Propdlg.m_nMs = m_nMs;         //Number of Inner Molecules
	if (m_nMs == 700) m_Propdlg.m_nMs = 0; 
	m_Propdlg.m_dt = int(m_fdt*1000);

	m_Propdlg.m_ndamping = m_ndamping;

	m_Propdlg.m_fMrad = int(m_fMrad*1000.0);         //Radius of Sphere
	m_Propdlg.m_fMdist = int(m_fMdist*1000.0);        //Distance between Molecules 
	//m_Propdlg.m_fH = int(m_fH*10.0);            //Width of the Channel
	m_Propdlg.m_nWtype = m_nWtype;            //Wall type 0: BCC  1:FCC
	m_Propdlg.m_ncutoff = m_ncutoff;
	m_Propdlg.m_nwallinteraction = m_nwallinteraction;
	m_Propdlg.m_nroop = m_nroop;
	m_Propdlg.m_nmilisec = m_nmilisec;
	m_Propdlg.m_targett = m_Ttarget; 
	m_Propdlg.m_tt = wallTempT;         
	m_Propdlg.m_bt = wallTempB;
	m_Propdlg.m_istrt=m_istrt;
	m_Propdlg.m_acc=m_acc;
	m_Propdlg.m_max=m_max;
	m_Propdlg.m_spring=m_spring;
	m_Propdlg.m_fH2=m_fH2;
	m_Propdlg.m_fW2=m_fW2;
	m_Propdlg.m_fZ2=m_fZ2;



	m_Propdlg.DoModal();

	m_istrt= m_Propdlg.m_istrt;
	m_dEwall = m_dEgas*m_istrt;

	m_acc= m_Propdlg.m_acc;
	
	m_spring= m_Propdlg.m_spring;
	d_K_M = m_spring*16.0*4.0*m_dEgas/(m_dSigma*m_dSigma);

	m_fH2= m_Propdlg.m_fH2;
	m_fH = 0.540*m_fH2*1.0;
	m_fW2= m_Propdlg.m_fW2;
	m_fW = 0.540*m_fW2*1.0;
	m_fZ2= m_Propdlg.m_fZ2;
	m_fZ = 0.540*m_fZ2*1.0;
	
	double tau;
	tau = m_dSigma*sqrt(4.0/m_dEgas); //= 2.1563481794417
	m_dGravP = m_acc*m_dSigma/(tau*tau);

	m_max= m_Propdlg.m_max;
	m_Umax = m_max*sqrt(m_dEgas/4.0);

	m_fdt = double(m_Propdlg.m_dt)/1000.0;

	m_nwMs = m_Propdlg.m_nwMs;           //Number of Wall Molecules
	m_nMs = m_Propdlg.m_nMs;         //Number of Inner Molecules

	m_fMrad = double(m_Propdlg.m_fMrad)*0.001;         //Radius of Sphere
//	m_fMdist = 	double(m_Propdlg.m_fMdist)*0.001;        //Distance between Molecules 
//	m_fH = double(m_Propdlg.m_fH)*0.1;            //Width of the Channel  Do Not modify here
	m_nWtype = m_Propdlg.m_nWtype;            //Wall type 0: BCC  1:FCC
	m_ncutoff = m_Propdlg.m_ncutoff;
	m_nwallinteraction = m_Propdlg.m_nwallinteraction;
	m_ndamping = m_Propdlg.m_ndamping;
	m_nroop = m_Propdlg.m_nroop;
	m_nmilisec = m_Propdlg.m_nmilisec;
	m_Ttarget = m_Propdlg.m_targett; 
	wallTempT = m_Propdlg.m_tt;         
	wallTempB = m_Propdlg.m_bt;

	//SetInitialCondition();	
	DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());

	
}

void CMFCGLView::OnAutoTimer()
{
	// TODO: Add your command handler code here

	omp_set_num_threads(iCPU);


	if (m_time == 0)
	{
		m_time =1;
		SetTimer(1,m_nmilisec,NULL);
		m_dautocolor = 0.8;
		m_Gcount = 0; //Set Chart
		DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
	}
	else
	{
		m_time=0;
		KillTimer(1);
		m_dautocolor = 0.4;
		DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
	}

}

void CMFCGLView::OnDamp1()
{
	// TODO: Add your command handler code here
	m_dampcoef = 1.0;

}

void CMFCGLView::OnDampdn()
{
		// TODO: Add your command handler code here
	m_dampcoef = m_dampcoef - 0.01;
}

void CMFCGLView::OnDampup()
{
		// TODO: Add your command handler code here
	m_dampcoef = m_dampcoef + 0.01;
}

void CMFCGLView::OnFileOpen()
{
	

/*	m_nMs = m_nwMs;

		for(ii = 0; ii<m_nMs ; ii++)//Top wall
	{
		if ( m_aM[ii].m_fX > (m_fH-0.3) || m_aM[ii].m_fY > (m_fH-0.3) )
		{
		m_aM[ii].m_fX = m_aBwM[ii].m_fX - m_fH+0.3 ;
		m_aM[ii].m_fY = m_aBwM[ii].m_fY + 2.0;
		m_aM[ii].m_fZ = m_aBwM[ii].m_fZ - m_fH +0.3;
		}else
		{
		m_aM[ii].m_fX = m_aBwM[ii].m_fX;
		m_aM[ii].m_fY = m_aBwM[ii].m_fY+1.0;
		m_aM[ii].m_fZ = m_aBwM[ii].m_fZ;
		}
	}
*/

//read fluid

	int ii;
	CFileDialog dlg(TRUE);
	dlg.DoModal();

	if (dlg.GetFileName() == "")
	{
		return;
	}

	m_rfile.Open(dlg.GetFileName(),CFile::modeRead);

	CString temp;
	m_rfile.ReadString(temp);
	
	m_nMs = atoi(temp);
//	m_nMs = 549; //Array starts from 0

	for( ii = 0; ii<m_nMs ; ii++) //Write x,y,vx,vy
	{

		
	char   *stopstring;

	m_rfile.ReadString(temp);
	m_aM[ii].m_fX = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);
	m_aM[ii].m_fY = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);
	m_aM[ii].m_fZ = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);

	m_aM[ii].m_fVx = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);
	m_aM[ii].m_fVy = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);
	m_aM[ii].m_fVz = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);

//	if (m_aM[ii].m_fY <= m_fH/3.0)
//	{m_aM[ii].m_fY = m_aM[ii].m_fY+0.15;}

//	if (m_aM[ii].m_fY >= 2*m_fH/3.0)
//	{m_aM[ii].m_fY = m_aM[ii].m_fY-0.15;}
///	m_aM[ii].m_fX = m_aM[ii].m_fX-0.15;
//	m_aM[ii].m_fY = m_aM[ii].m_fY-0.15;
//	m_aM[ii].m_fZ = m_aM[ii].m_fZ-0.15;
	
	}

	m_rfile.Close();
	
SortGridInit();



//read wall cancel 

	/*
//	CFileDialog dlg(TRUE);
	dlg.DoModal();

	if (dlg.GetFileName() == "fluid")
	{
		return;
	}

	ifstream wallfilein(dlg.GetFileName());

	wallfilein >> m_nwMs;

	int dump;

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		wallfilein >> dump;
		wallfilein >> m_aBwM[ii].m_fX;
		wallfilein >> m_aBwM[ii].m_fY;
		wallfilein >> m_aBwM[ii].m_fZ;

		m_aBwM[ii].m_fY = m_aBwM[ii].m_fY - m_fMdist; //for m_fh clearance.
		//m_aBwM[ii].m_fX = m_aBwM[ii].m_fX - m_fMdist; 
		//m_aBwM[ii].m_fZ = m_aBwM[ii].m_fZ - m_fMdist; 

	}
		
	wallfilein.close();

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		
		m_aTwM[ii].m_fX = m_aBwM[ii].m_fX;//-m_fMdist;
		m_aTwM[ii].m_fY = m_aBwM[ii].m_fY + m_fH + m_fMdist;
		m_aTwM[ii].m_fZ = m_aBwM[ii].m_fZ;//-m_fMdist;

	}
	*/

//	OnRun();
}

void CMFCGLView::OnFileSave()
{
	
	CFileDialog dlg(FALSE);
	dlg.DoModal();

	if (dlg.GetFileName() == "")
	{
		return;
	}


	m_file.Open(dlg.GetFileName(),CFile::modeCreate | CFile::modeWrite);

	CString temp;
//	char buffer[20];

	//_itoa(m_nMs, buffer, 10 );
	//temp = CString(buffer);
	temp.Format("%d",m_nMs);
	m_file.WriteString(temp+"\n");  //Write Number of Moles

//	char buffer1[20],buffer2[20],buffer3[20],buffer4[20],buffer5[20],buffer6[20];

	for(int ii = 0; ii<m_nMs ; ii++) //Write x,y,vx,vy
	{
				

				//_gcvt( m_aM[ii].m_fX, 18, buffer1);
				//m_file.WriteString(buffer1);
				temp.Format("%lf",m_aM[ii].m_fX);
				m_file.WriteString(temp);
				m_file.WriteString("\n"); 
	
				//_gcvt( m_aM[ii].m_fY, 18, buffer2);
				//m_file.WriteString(buffer2); 
				temp.Format("%lf",m_aM[ii].m_fY);
				m_file.WriteString(temp);
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fZ, 18, buffer5);
				//m_file.WriteString(buffer5); 
				temp.Format("%lf",m_aM[ii].m_fZ);
				m_file.WriteString(temp);
				m_file.WriteString("\n"); 


				//_gcvt( m_aM[ii].m_fVx, 18, buffer3);
				//m_file.WriteString(buffer3);
				temp.Format("%lf",m_aM[ii].m_fVx);
				m_file.WriteString(temp);
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fVy, 18, buffer4);
				//m_file.WriteString(buffer4); 
				temp.Format("%lf",m_aM[ii].m_fVy);
				m_file.WriteString(temp);
				m_file.WriteString("\n"); 

				//_gcvt( m_aM[ii].m_fVz, 18, buffer6);
				//m_file.WriteString(buffer6); 
				temp.Format("%lf",m_aM[ii].m_fVz);
				m_file.WriteString(temp);
				m_file.WriteString("\n"); 

				m_file.WriteString("\n"); 
	}


	m_file.Close();
	
	  
}

void CMFCGLView::OnWalldown()
{
		if (m_nWallDown == 0)	
	{
		m_nWallDown = 1;
		d_vibw = 1.0; 
	}else
	{
		m_nWallDown = 0;
		d_vibw = 0.0;
	}
 	
}

void CMFCGLView::OnReverseflot()
{
	if (m_rf == 1.0) m_rf = -1.0;
	else m_rf = 1.0;

	//DrawScene(); SwapBuffers(m_pDC->GetSafeHdc());
}

void CMFCGLView::OnBoltzmannVel()
{
	
	double V1,V2,V3,SD;
	double A1,A2,A3;
	int ii;


		for(ii = 0; ii<m_nMs ; ii++)
		{
				//m_Ttarget = 120;
				SD = sqrt((1.0/m_Kb)*m_Ttarget);

				A1 = 0.0; A2 = 0.0; A3 = 0.0;

				for (int rr=0;rr<12;rr++)
				{
					A1 = A1 + double(rand())/double(RAND_MAX);
					A2 = A2 + double(rand())/double(RAND_MAX);
					A3 = A3 + double(rand())/double(RAND_MAX);
				}
				V1 = ( A1 - 6.0 ) * SD;
				V2 = ( A2 - 6.0 ) * SD;
				V3 = ( A3 - 6.0 ) * SD;


				m_aM[ii].m_fVx = -V1;
				m_aM[ii].m_fVy = -V2;
				m_aM[ii].m_fVz = -V3;
		}

//for zero thermal velocity

		double sumx,sumy,sumz;
		sumx = 0;
		sumy = 0;
		sumz = 0;

		INT IsVelNonZero;
		IsVelNonZero = 1;

		while (IsVelNonZero)
		{
			for(ii = 0; ii<m_nMs ; ii++)
			{
					sumx = sumx + m_aM[ii].m_fVx;
					sumy = sumy + m_aM[ii].m_fVy;
					sumz = sumz + m_aM[ii].m_fVz;
			}


				if (fabs(sumx) > 1E-15)
				{
					for(ii = 0; ii<m_nMs ; ii++)
					{
					m_aM[ii].m_fVx = m_aM[ii].m_fVx - sumx/double(m_nMs);
					}
				}
				if (fabs(sumy) > 1E-15)
				{
					for(ii = 0; ii<m_nMs ; ii++)
					{
					m_aM[ii].m_fVy = m_aM[ii].m_fVy - sumy/double(m_nMs);
					}
				}
				if (fabs(sumz) > 1E-15)
				{
					for(ii = 0; ii<m_nMs ; ii++)
					{
					m_aM[ii].m_fVz = m_aM[ii].m_fVz - sumz/double(m_nMs);
					}
				}


			for( ii = 0; ii<m_nMs ; ii++)
			{
				m_V = m_V +  pow(m_aM[ii].m_fVx,2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
			} 

			m_T = (1.0/3.0)*m_Kb*m_V/double(m_nMs);
			m_V = 0;

				for(ii = 0; ii<m_nMs ; ii++)
				{
						m_aM[ii].m_fVx = m_aM[ii].m_fVx *sqrt(m_Ttarget/m_T);
						m_aM[ii].m_fVy = m_aM[ii].m_fVy *sqrt(m_Ttarget/m_T) ;
						m_aM[ii].m_fVz = m_aM[ii].m_fVz *sqrt(m_Ttarget/m_T) ;
				}

			for( ii = 0; ii<m_nMs ; ii++)
			{
				m_V = m_V +  pow(m_aM[ii].m_fVx,2.0) + pow(m_aM[ii].m_fVy,2.0) + pow(m_aM[ii].m_fVz,2.0);
			} 

			m_T = (1.0/3.0)*m_Kb*m_V/double(m_nMs);
			m_V = 0;
	
				sumx = 0;
				sumy = 0;
				sumz = 0;


			for(ii = 0; ii<m_nMs ; ii++)
			{
					sumx = sumx + m_aM[ii].m_fVx;
					sumy = sumy + m_aM[ii].m_fVy;
					sumz = sumz + m_aM[ii].m_fVz;
			}

				
		double difftemp;

			difftemp = m_Ttarget-m_T;
				if (fabs(sumx) < 1E-15 && fabs(sumy) < 1E-15 && fabs(sumz) < 1E-15)
				{
					IsVelNonZero = 0;
				}
				if ( fabs(difftemp) > 1E-13) 
				{
					IsVelNonZero = 1;
				}
				sumx = 0;
				sumy = 0;
				sumz = 0;

		}

//for zero thermal velocity


	//Boltzmann for B wall

		

		for(ii = 0; ii<m_nwMs ; ii++)
		{
				//m_Ttarget = 120;
				SD = sqrt((1.0/m_Kb)*m_Ttarget);

				A1 = 0.0; A2 = 0.0; A3 = 0.0;

				for (int rr=0;rr<12;rr++)
				{
					A1 = A1 + double(rand())/double(RAND_MAX);
					A2 = A2 + double(rand())/double(RAND_MAX);
					A3 = A3 + double(rand())/double(RAND_MAX);
				}
				V1 = ( A1 - 6.0 ) * SD;
				V2 = ( A2 - 6.0 ) * SD;
				V3 = ( A3 - 6.0 ) * SD;


				m_aBwM[ii].m_fVx = -V1;
				m_aBwM[ii].m_fVy = -V2;
				m_aBwM[ii].m_fVz = -V3;
		}


		
		sumx = 0;
		sumy = 0;
		sumz = 0;
		IsVelNonZero = 1;

		while (IsVelNonZero)
		{
			for(ii = 0; ii<m_nwMs ; ii++)
			{
					sumx = sumx + m_aBwM[ii].m_fVx;
					sumy = sumy + m_aBwM[ii].m_fVy;
					sumz = sumz + m_aBwM[ii].m_fVz;
			}


				if (fabs(sumx) > 1E-15)
				{
					for(ii = 0; ii<m_nwMs ; ii++)
					{
					m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx - sumx/double(m_nwMs);
					}
				}
				if (fabs(sumy) > 1E-15)
				{
					for(ii = 0; ii<m_nwMs ; ii++)
					{
					m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy - sumy/double(m_nwMs);
					}
				}
				if (fabs(sumz) > 1E-15)
				{
					for(ii = 0; ii<m_nwMs ; ii++)
					{
					m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz - sumz/double(m_nwMs);
					}
				}


			for( ii = 0; ii<m_nwMs ; ii++)
			{
				m_V = m_V +  pow(m_aBwM[ii].m_fVx,2.0) + pow(m_aBwM[ii].m_fVy,2.0) + pow(m_aBwM[ii].m_fVz,2.0);
			} 

			m_T = (1.0/3.0)*m_Kb*m_V/double(m_nwMs);
			m_V = 0;

				for(ii = 0; ii<m_nwMs ; ii++)
				{
						m_aBwM[ii].m_fVx = m_aBwM[ii].m_fVx *sqrt(m_Ttarget/m_T);
						m_aBwM[ii].m_fVy = m_aBwM[ii].m_fVy *sqrt(m_Ttarget/m_T) ;
						m_aBwM[ii].m_fVz = m_aBwM[ii].m_fVz *sqrt(m_Ttarget/m_T) ;
				}

			for( ii = 0; ii<m_nwMs ; ii++)
			{
				m_V = m_V +  pow(m_aBwM[ii].m_fVx,2.0) + pow(m_aBwM[ii].m_fVy,2.0) + pow(m_aBwM[ii].m_fVz,2.0);
			} 

			m_T = (1.0/3.0)*m_Kb*m_V/double(m_nwMs);
			m_V = 0;
	
				sumx = 0;
				sumy = 0;
				sumz = 0;


			for(ii = 0; ii<m_nwMs ; ii++)
			{
					sumx = sumx + m_aBwM[ii].m_fVx;
					sumy = sumy + m_aBwM[ii].m_fVy;
					sumz = sumz + m_aBwM[ii].m_fVz;
			}

				
		double difftemp;

			difftemp = m_Ttarget-m_T;
				if (fabs(sumx) < 1E-15 && fabs(sumy) < 1E-15 && fabs(sumz) < 1E-15)
				{
					IsVelNonZero = 0;
				}
				if ( fabs(difftemp) > 1E-13) 
				{
					IsVelNonZero = 1;
				}
				sumx = 0;
				sumy = 0;
				sumz = 0;

		}


	//Boltzmann for T wall

		

		for(ii = 0; ii<m_nwMs ; ii++)
		{
				//m_Ttarget = 120;
				SD = sqrt((1.0/m_Kb)*m_Ttarget);

				A1 = 0.0; A2 = 0.0; A3 = 0.0;

				for (int rr=0;rr<12;rr++)
				{
					A1 = A1 + double(rand())/double(RAND_MAX);
					A2 = A2 + double(rand())/double(RAND_MAX);
					A3 = A3 + double(rand())/double(RAND_MAX);
				}
				V1 = ( A1 - 6.0 ) * SD;
				V2 = ( A2 - 6.0 ) * SD;
				V3 = ( A3 - 6.0 ) * SD;


				m_aTwM[ii].m_fVx = -V1;
				m_aTwM[ii].m_fVy = -V2;
				m_aTwM[ii].m_fVz = -V3;
		}


		
		sumx = 0;
		sumy = 0;
		sumz = 0;
		IsVelNonZero = 1;

		while (IsVelNonZero)
		{
			for(ii = 0; ii<m_nwMs ; ii++)
			{
					sumx = sumx + m_aTwM[ii].m_fVx;
					sumy = sumy + m_aTwM[ii].m_fVy;
					sumz = sumz + m_aTwM[ii].m_fVz;
			}


				if (fabs(sumx) > 1E-15)
				{
					for(ii = 0; ii<m_nwMs ; ii++)
					{
					m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx - sumx/double(m_nwMs);
					}
				}
				if (fabs(sumy) > 1E-15)
				{
					for(ii = 0; ii<m_nwMs ; ii++)
					{
					m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy - sumy/double(m_nwMs);
					}
				}
				if (fabs(sumz) > 1E-15)
				{
					for(ii = 0; ii<m_nwMs ; ii++)
					{
					m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz - sumz/double(m_nwMs);
					}
				}


			for( ii = 0; ii<m_nwMs ; ii++)
			{
				m_V = m_V +  pow(m_aTwM[ii].m_fVx,2.0) + pow(m_aTwM[ii].m_fVy,2.0) + pow(m_aTwM[ii].m_fVz,2.0);
			} 

			m_T = (1.0/3.0)*m_Kb*m_V/double(m_nwMs);
			m_V = 0;

				for(ii = 0; ii<m_nwMs ; ii++)
				{
						m_aTwM[ii].m_fVx = m_aTwM[ii].m_fVx *sqrt(m_Ttarget/m_T);
						m_aTwM[ii].m_fVy = m_aTwM[ii].m_fVy *sqrt(m_Ttarget/m_T) ;
						m_aTwM[ii].m_fVz = m_aTwM[ii].m_fVz *sqrt(m_Ttarget/m_T) ;
				}

			for( ii = 0; ii<m_nwMs ; ii++)
			{
				m_V = m_V +  pow(m_aTwM[ii].m_fVx,2.0) + pow(m_aTwM[ii].m_fVy,2.0) + pow(m_aTwM[ii].m_fVz,2.0);
			} 

			m_T = (1.0/3.0)*m_Kb*m_V/double(m_nwMs);
			m_V = 0;
	
				sumx = 0;
				sumy = 0;
				sumz = 0;


			for(ii = 0; ii<m_nwMs ; ii++)
			{
					sumx = sumx + m_aTwM[ii].m_fVx;
					sumy = sumy + m_aTwM[ii].m_fVy;
					sumz = sumz + m_aTwM[ii].m_fVz;
			}

				
		double difftemp;

			difftemp = m_Ttarget-m_T;
				if (fabs(sumx) < 1E-15 && fabs(sumy) < 1E-15 && fabs(sumz) < 1E-15)
				{
					IsVelNonZero = 0;
				}
				if ( fabs(difftemp) > 1E-13) 
				{
					IsVelNonZero = 1;
				}
				sumx = 0;
				sumy = 0;
				sumz = 0;

		}



}



void CMFCGLView::OnOpenwall()
{
		int ii;
	CFileDialog dlg(TRUE);
	dlg.DoModal();

	if (dlg.GetFileName() == "")
	{
		return;
	}


	m_rfile.Open(dlg.GetFileName(),CFile::modeRead);
	
	CString temp;
	m_rfile.ReadString(temp);

	m_nwMs = atoi(temp);
	//wallfilein >> m_nwMs;
	m_rfile.ReadString(temp); //1st blank

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{

		char   *stopstring;

		m_rfile.ReadString(temp); //numbering	

		m_rfile.ReadString(temp);
		m_aBwM[ii].m_fX = strtod(temp,&stopstring);
		m_rfile.ReadString(temp);
		m_aBwM[ii].m_fY = strtod(temp,&stopstring);
		m_rfile.ReadString(temp);
		m_aBwM[ii].m_fZ = strtod(temp,&stopstring);
		m_rfile.ReadString(temp); //end blank

		m_aBwM[ii].m_fY = m_aBwM[ii].m_fY - m_fMdist; //for m_fh clearance.
		//m_aBwM[ii].m_fX = m_aBwM[ii].m_fX - m_fMdist; 
		//m_aBwM[ii].m_fZ = m_aBwM[ii].m_fZ - m_fMdist; 

	}


	m_rfile.Close();

/* this wiht using namespace std might cause error in vs2008
	ifstream wallfilein(dlg.GetFileName());

	wallfilein >> m_nwMs;

	int dump;

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		wallfilein >> dump;
		wallfilein >> m_aBwM[ii].m_fX;
		wallfilein >> m_aBwM[ii].m_fY;
		wallfilein >> m_aBwM[ii].m_fZ;

		m_aBwM[ii].m_fY = m_aBwM[ii].m_fY - m_fMdist; //for m_fh clearance.
		//m_aBwM[ii].m_fX = m_aBwM[ii].m_fX - m_fMdist; 
		//m_aBwM[ii].m_fZ = m_aBwM[ii].m_fZ - m_fMdist; 

	}
		
	wallfilein.close();


*/



	//remove the boundary molecules// m_fW,FZ

	int tmpcount;
	tmpcount = 0;
	
	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		
		m_aBwM[tmpcount].m_fX = m_aBwM[ii].m_fX;
		m_aBwM[tmpcount].m_fY = m_aBwM[ii].m_fY;
		m_aBwM[tmpcount].m_fZ = m_aBwM[ii].m_fZ;
		if(m_aBwM[ii].m_fX >= m_fW)	{tmpcount = tmpcount - 1;}
		tmpcount = tmpcount +1;
	}

	m_nwMs = tmpcount;

	tmpcount = 0;
	
	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		m_aBwM[tmpcount].m_fX = m_aBwM[ii].m_fX;
		m_aBwM[tmpcount].m_fY = m_aBwM[ii].m_fY;
		m_aBwM[tmpcount].m_fZ = m_aBwM[ii].m_fZ;
		if(m_aBwM[ii].m_fZ >= m_fZ)	{tmpcount = tmpcount - 1;}
		tmpcount = tmpcount +1;
	}

	m_nwMs = tmpcount;

	//remove the boundary molecules// m_fW,FZ


	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		
		m_aTwM[ii].m_fX = m_aBwM[ii].m_fX;//-m_fMdist;
		m_aTwM[ii].m_fY = m_aBwM[ii].m_fY + m_fH + m_fMdist;
		m_aTwM[ii].m_fZ = m_aBwM[ii].m_fZ;//-m_fMdist;

	}


	////////////////////////Initialize wall molecules

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{

		m_aBwM[ii].m_fXo = m_aBwM[ii].m_fX;
		m_aBwM[ii].m_fYo = m_aBwM[ii].m_fY;
		m_aBwM[ii].m_fZo = m_aBwM[ii].m_fZ;
	
		m_aTwM[ii].m_fXo = m_aTwM[ii].m_fX;
		m_aTwM[ii].m_fYo = m_aTwM[ii].m_fY;
		m_aTwM[ii].m_fZo = m_aTwM[ii].m_fZ;

		m_aBwM[ii].m_fVx = 0.0;
		m_aBwM[ii].m_fVy = 0.0;
		m_aBwM[ii].m_fVz = 0.0;
	
		m_aTwM[ii].m_fVx = 0.0;
		m_aTwM[ii].m_fVy = 0.0;
		m_aTwM[ii].m_fVz = 0.0;


	
	}



	SortGridInit();
//	SortWallMolecules();
	//OnRun();
}

void CMFCGLView::OnMsave()
{
		int Multiplex,Multipley,Multiplez,m_nMsTemp;
	double fx,fy,fz;
	Multiplex = 1;
	Multipley = 4;
	Multiplez = 1;


	CFileDialog dlg(FALSE);
	dlg.DoModal();

	if (dlg.GetFileName() == "")
	{
		return;
	}

	m_nMsTemp = m_nMs*Multiplex*Multipley*Multiplez;


	m_file.Open(dlg.GetFileName(),CFile::modeCreate | CFile::modeWrite);

	CString temp;
	//char buffer[20];

	//_itoa(m_nMsTemp, buffer, 10 );
	//temp = CString(buffer);

	temp.Format("%d",m_nMsTemp);

	m_file.WriteString(temp+"\n");  //Write Number of Moles

	//char buffer1[20],buffer2[20],buffer3[20],buffer4[20],buffer5[20],buffer6[20];

	for (int xx=0; xx < Multiplex ; xx++)
	{
			for (int yy=0; yy < Multipley ; yy++)
		{
					for (int zz=0; zz < Multiplez ; zz++)
			{

						for(int ii = 0; ii<m_nMs ; ii++) //Write x,y,vx,vy
						{
							fx = m_aM[ii].m_fX + double(xx)*m_fW;
							//_gcvt( fx, 18, buffer1);
							//m_file.WriteString(buffer1); 
							temp.Format("%lf",fx);	m_file.WriteString(temp);
							m_file.WriteString("\n"); 

							fy = m_aM[ii].m_fY + double(yy)*m_fH;
							//_gcvt( fy, 18, buffer2);
							//m_file.WriteString(buffer2); 
							temp.Format("%lf",fy);	m_file.WriteString(temp);
							m_file.WriteString("\n"); 
							
							fz = m_aM[ii].m_fZ + double(zz)*m_fZ;
							//_gcvt( fz, 18, buffer5);
							//m_file.WriteString(buffer5); 
							temp.Format("%lf",fz);	m_file.WriteString(temp);
							m_file.WriteString("\n"); 

							//_gcvt( m_aM[ii].m_fVx, 18, buffer3);
							//m_file.WriteString(buffer3); 
							temp.Format("%lf",m_aM[ii].m_fVx);	m_file.WriteString(temp);
							m_file.WriteString("\n"); 

							//_gcvt( m_aM[ii].m_fVy, 18, buffer4);
							//m_file.WriteString(buffer4); 
							temp.Format("%lf",m_aM[ii].m_fVy);	m_file.WriteString(temp);
							m_file.WriteString("\n"); 

							//_gcvt( m_aM[ii].m_fVz, 18, buffer6);
							//m_file.WriteString(buffer6); 
							temp.Format("%lf",m_aM[ii].m_fVz);	m_file.WriteString(temp);
							m_file.WriteString("\n"); 
							
							
							m_file.WriteString("\n"); 
						}//
			}//zz
		}//yy
	}//for (int xx=0; xx < Multiple < xx+1)



	m_file.Close();
	
}

void CMFCGLView::OnPost()
{
	


		double ramda,N_totstep,M_avgstep;
		N_totstep = 100000;
		M_avgstep = 8000;
		int N,M;
		ramda = 0.0;
		N=0; M=0;


	int ii;
	CFileDialog dlg(TRUE);
	dlg.DoModal();

	if (dlg.GetFileName() == "")
	{
		return;
	}

	m_rfile.Open(dlg.GetFileName(),CFile::modeRead);

	CString temp;
	m_rfile.ReadString(temp);
	
	m_nMs = atoi(temp);
//	m_nMs = 549; //Array starts from 0

	for( ii = 0; ii<N_totstep ; ii++) //Write x,y,vx,vy
	{
	char   *stopstring;
	m_rfile.ReadString(temp);
	d_J[ii][0] = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);
	d_J[ii][1] = strtod(temp,&stopstring);
	m_rfile.ReadString(temp);
	d_J[ii][2] = strtod(temp,&stopstring);
	}

	m_rfile.Close();
	//double d_J[100000][3];
	//int  	n_Jcount,n_Jtotal;

		for( M = 1; M < M_avgstep+1 ; M++)
		{
			for( N = 1; N < N_totstep-M_avgstep+1 ; N++)
			{
			ramda = ramda + d_J[M+N-1][0]*d_J[N-1][0] + d_J[M+N-1][1]*d_J[N-1][1] + d_J[M+N-1][2]*d_J[N-1][2]; 
			}
		}

		ramda =ramda*(m_fdt/s_to_ps)/(3.0*((m_fH*m_fW*m_fZ)/(m_to_nm*m_to_nm*m_to_nm))*(m_BoltzN/joule_to_gnmps)*m_Ttarget*m_Ttarget*(N_totstep-M_avgstep));
		ramda = 0;
		//ramda should be around 0.11 w/mk
}


BOOL CMFCGLView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	// TODO: 여기에 메시지 처리기 코드를 추가 및/또는 기본값을 호출합니다.

	 if (zDelta > 0)
	 {
		m_zTranslate=m_zTranslate + 1.0;
	 }else{
		m_zTranslate=m_zTranslate - 1.0;
	 }
	DrawScene();
	SwapBuffers(m_pDC->GetSafeHdc());



	return CView::OnMouseWheel(nFlags, zDelta, pt);
}

void CMFCGLView::visualclose()
{	if (visual != 1) visual = 1;
	else visual = 0;
	return;
}
void CMFCGLView::OnLattice()
{
	m_LatticeDlg.m_nx = m_nx;
	m_LatticeDlg.m_ny = m_ny;
	m_LatticeDlg.m_nz = m_nz;

	m_LatticeDlg.DoModal();

	m_nx = m_LatticeDlg.m_nx;
	m_ny = m_LatticeDlg.m_ny;
	m_nz = m_LatticeDlg.m_nz;

		double Hcoef;
	Hcoef = m_fMdist*(2.0/sqrt(2.0)); //0.54022958082652
	Hcoef = 0.540;
	//double offsetx;
	//double offsety;
	//double offsetz;
//m_dSigma = 0.3405; 
//m_fH = 0.540*6.0*1.0;   

	//left side wall
	m_aM[0].m_fX = 0.0;
	m_aM[0].m_fY = 0.0;
	m_aM[0].m_fZ = 0.0;

	m_aM[1].m_fX = 0.0;
	m_aM[1].m_fY = Hcoef;
	m_aM[1].m_fZ = 0.0;

	m_aM[2].m_fX = 0.0;
	m_aM[2].m_fY = Hcoef/2;
	m_aM[2].m_fZ = Hcoef/2;

//=========================================================

	//front side wall
	m_aM[3].m_fX = 0.0;
	m_aM[3].m_fY = 0.0;
	m_aM[3].m_fZ = Hcoef;

	m_aM[4].m_fX = 0.0;
	m_aM[4].m_fY = Hcoef;
	m_aM[4].m_fZ = Hcoef;

//=========================================================

	m_aM[5].m_fX = Hcoef;
	m_aM[5].m_fY = 0.0;
	m_aM[5].m_fZ = Hcoef;

	m_aM[6].m_fX = Hcoef;
	m_aM[6].m_fY = Hcoef;
	m_aM[6].m_fZ = Hcoef;

	m_aM[7].m_fX = Hcoef/2;
	m_aM[7].m_fY = Hcoef/2;
	m_aM[7].m_fZ = Hcoef;

//===========================================================

	m_aM[8].m_fX = Hcoef;
	m_aM[8].m_fY = 0.0;
	m_aM[8].m_fZ = 0.0;

	m_aM[9].m_fX = Hcoef;
	m_aM[9].m_fY = Hcoef;
	m_aM[9].m_fZ = 0.0;

	m_aM[10].m_fX = Hcoef/2;
	m_aM[10].m_fY = 0.0;
	m_aM[10].m_fZ = Hcoef/2;

	m_aM[11].m_fX = Hcoef/2;
	m_aM[11].m_fY = Hcoef/2;
	m_aM[11].m_fZ = 0.0;

	m_aM[12].m_fX = Hcoef;
	m_aM[12].m_fY = Hcoef/2;
	m_aM[12].m_fZ = Hcoef/2;

	m_aM[13].m_fX = Hcoef/2;
	m_aM[13].m_fY = Hcoef;
	m_aM[13].m_fZ = Hcoef/2;

	int n_x;
	int n_z;
	int n_y;
	int xxi,zzi,yyi,rptn;

	n_x = m_nx;
	n_z = m_nz;
	n_y = m_ny;


	n_x = n_x - 1;
	n_z = n_z - 1;
	n_y = n_y - 1;
	rptn = 14;

	for (xxi=0; xxi < n_x ; xxi++)
	{
		for (int n = 0 ; n < 9 ; n++)
		{
			m_aM[rptn].m_fX = m_aM[n+5].m_fX + Hcoef*(xxi+1);
			m_aM[rptn].m_fY = m_aM[n+5].m_fY;
			m_aM[rptn].m_fZ = m_aM[n+5].m_fZ;

			rptn = rptn +1;
		}
	}


	for (zzi=0; zzi < n_z ; zzi++)
	{
		for (int n = 0 ; n < (5 + 9*(n_x+1))  ; n++)
		{
			if (m_aM[n].m_fZ != 0)
			{
			m_aM[rptn].m_fX = m_aM[n].m_fX;
			m_aM[rptn].m_fY = m_aM[n].m_fY;
			m_aM[rptn].m_fZ = m_aM[n].m_fZ + Hcoef*(zzi+1);
			rptn = rptn +1;
			}
		}
	}
	
	int tem_number;
	tem_number = rptn;


	for (yyi=0; yyi < n_y ; yyi++)
	{
		for (int n = 0 ; n < tem_number  ; n++)
		{
			if (m_aM[n].m_fY != 0)
			{
			m_aM[rptn].m_fX = m_aM[n].m_fX;
			m_aM[rptn].m_fY = m_aM[n].m_fY + Hcoef*(yyi+1);
			m_aM[rptn].m_fZ = m_aM[n].m_fZ;
			rptn = rptn +1;
			}
		}
	}


	//let's erase duplicated molecules


	tem_number = rptn;


		for (int ii = 0 ; ii < rptn  ; ii++)
		{
			for (int jj = 0 ; jj < rptn  ; jj++)
			{
				if (ii!=jj)
				{
					if (m_aM[ii].m_fX == m_aM[jj].m_fX)
					{
						if (m_aM[ii].m_fY == m_aM[jj].m_fY)
						{
							if (m_aM[ii].m_fZ == m_aM[jj].m_fZ)
							{
								AfxMessageBox("duplicate");
						//		for (int kk = jj ; kk < tem_number  ; kk++)
						//		{
						//			m_aBwM[kk].m_fX = m_aBwM[kk+1].m_fX;
							//		m_aBwM[kk].m_fY = m_aBwM[kk+1].m_fY;
							//		m_aBwM[kk].m_fZ = m_aBwM[kk+1].m_fZ;
							//	}

							//	rptn = rptn - 1;
							}
						}
					}
				}

			}//for (int jj = 0 ; jj < tem_number  ; jj++)
		}





	m_nMs = rptn; //because the numbering of array starts from 0
	

	for(int ii = 0; ii<m_nMs ; ii++)//Top wall
	{
		m_aM[ii].m_fY = m_aM[ii].m_fY + m_dSigma;
	}



	OnRun();
			
	DrawScene();
	SwapBuffers(m_pDC->GetSafeHdc());

}

void CMFCGLView::OnS()
{

//////////////////////////////////BW
	m_SDlg.m_nnx = m_nnx;
	m_SDlg.m_nnz = m_nnz;

		m_SDlg.DoModal();

	m_nnx = m_SDlg.m_nnx;
	m_nnz = m_SDlg.m_nnz;

		double Hcoef;
	Hcoef = m_fMdist*(2.0/sqrt(2.0)); //0.54022958082652
	Hcoef = 0.540;
	//double offsetx;
	//double offsety;
	//double offsetz;


	//left side wall
	m_aBwM[0].m_fX = 0.0;
	m_aBwM[0].m_fY = 0.0-Hcoef;
	m_aBwM[0].m_fZ = 0.0;

	m_aBwM[1].m_fX = 0.0;
	m_aBwM[1].m_fY = Hcoef-Hcoef;
	m_aBwM[1].m_fZ = 0.0;

	m_aBwM[2].m_fX = 0.0;
	m_aBwM[2].m_fY = Hcoef/2-Hcoef;
	m_aBwM[2].m_fZ = Hcoef/2;

//=========================================================

	//front side wall
	m_aBwM[3].m_fX = 0.0;
	m_aBwM[3].m_fY = 0.0-Hcoef;
	m_aBwM[3].m_fZ = Hcoef;

	m_aBwM[4].m_fX = 0.0;
	m_aBwM[4].m_fY = Hcoef-Hcoef;
	m_aBwM[4].m_fZ = Hcoef;

//=========================================================

	m_aBwM[5].m_fX = Hcoef;
	m_aBwM[5].m_fY = 0.0-Hcoef;
	m_aBwM[5].m_fZ = Hcoef;

	m_aBwM[6].m_fX = Hcoef;
	m_aBwM[6].m_fY = Hcoef-Hcoef;
	m_aBwM[6].m_fZ = Hcoef;

	m_aBwM[7].m_fX = Hcoef/2;
	m_aBwM[7].m_fY = Hcoef/2-Hcoef;
	m_aBwM[7].m_fZ = Hcoef;

//===========================================================

	m_aBwM[8].m_fX = Hcoef;
	m_aBwM[8].m_fY = 0.0-Hcoef;
	m_aBwM[8].m_fZ = 0.0;

	m_aBwM[9].m_fX = Hcoef;
	m_aBwM[9].m_fY = Hcoef-Hcoef;
	m_aBwM[9].m_fZ = 0.0;

	m_aBwM[10].m_fX = Hcoef/2;
	m_aBwM[10].m_fY = 0.0-Hcoef;
	m_aBwM[10].m_fZ = Hcoef/2;

	m_aBwM[11].m_fX = Hcoef/2;
	m_aBwM[11].m_fY = Hcoef/2-Hcoef;
	m_aBwM[11].m_fZ = 0.0;

	m_aBwM[12].m_fX = Hcoef;
	m_aBwM[12].m_fY = Hcoef/2-Hcoef;
	m_aBwM[12].m_fZ = Hcoef/2;

	m_aBwM[13].m_fX = Hcoef/2;
	m_aBwM[13].m_fY = Hcoef-Hcoef;
	m_aBwM[13].m_fZ = Hcoef/2;

	int n_x;
	int n_z;
	int n_y;
	int xxi,zzi,yyi,rptn;

	n_x = m_nnx;
	n_z = m_nnz;
	n_y = m_nny;


	n_x = n_x - 1;
	n_z = n_z - 1;
	n_y = n_y - 1;
	rptn = 14;

	for (xxi=0; xxi < n_x ; xxi++)
	{
		for (int n = 0 ; n < 9 ; n++)
		{
			m_aBwM[rptn].m_fX = m_aBwM[n+5].m_fX + Hcoef*(xxi+1);
			m_aBwM[rptn].m_fY = m_aBwM[n+5].m_fY;
			m_aBwM[rptn].m_fZ = m_aBwM[n+5].m_fZ;

			rptn = rptn +1;
		}
	}


	for (zzi=0; zzi < n_z ; zzi++)
	{
		for (int n = 0 ; n < (5 + 9*(n_x+1))  ; n++)
		{
			if (m_aBwM[n].m_fZ != 0)
			{
			m_aBwM[rptn].m_fX = m_aBwM[n].m_fX;
			m_aBwM[rptn].m_fY = m_aBwM[n].m_fY;
			m_aBwM[rptn].m_fZ = m_aBwM[n].m_fZ + Hcoef*(zzi+1);
			rptn = rptn +1;
			}
		}
	}
	
	int tem_number;
	tem_number = rptn;


	for (yyi=0; yyi < n_y ; yyi++)
	{
		for (int n = 0 ; n < tem_number  ; n++)
		{
			if (m_aBwM[n].m_fY != 0)
			{
			m_aBwM[rptn].m_fX = m_aBwM[n].m_fX;
			m_aBwM[rptn].m_fY = m_aBwM[n].m_fY + Hcoef*(yyi+1);
			m_aBwM[rptn].m_fZ = m_aBwM[n].m_fZ;
			rptn = rptn +1;
			}
		}
	}


	//let's erase duplicated molecules


	tem_number = rptn;


		for (int ii = 0 ; ii < rptn  ; ii++)
		{
			for (int jj = 0 ; jj < rptn  ; jj++)
			{
				if (ii!=jj)
				{
					if (m_aBwM[ii].m_fX == m_aBwM[jj].m_fX)
					{
						if (m_aBwM[ii].m_fY == m_aBwM[jj].m_fY)
						{
							if (m_aBwM[ii].m_fZ == m_aBwM[jj].m_fZ)
							{
								AfxMessageBox("duplicate");
						//		for (int kk = jj ; kk < tem_number  ; kk++)
						//		{
						//			m_aBwM[kk].m_fX = m_aBwM[kk+1].m_fX;
							//		m_aBwM[kk].m_fY = m_aBwM[kk+1].m_fY;
							//		m_aBwM[kk].m_fZ = m_aBwM[kk+1].m_fZ;
							//	}

							//	rptn = rptn - 1;
							}
						}
					}
				}

			}//for (int jj = 0 ; jj < tem_number  ; jj++)
		}




	m_nwMs = rptn; //because the numbering of array starts from 0
	

	int ii;

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{
		
		m_aTwM[ii].m_fX = m_aBwM[ii].m_fX;//-m_fMdist;
		m_aTwM[ii].m_fY = m_aBwM[ii].m_fY + m_fH + m_fMdist;
		m_aTwM[ii].m_fZ = m_aBwM[ii].m_fZ;//-m_fMdist;

	}


	////////////////////////Initialize wall molecules

	for(ii = 0; ii<m_nwMs ; ii++)//Top wall
	{

		m_aBwM[ii].m_fXo = m_aBwM[ii].m_fX;
		m_aBwM[ii].m_fYo = m_aBwM[ii].m_fY;
		m_aBwM[ii].m_fZo = m_aBwM[ii].m_fZ;
	
		m_aTwM[ii].m_fXo = m_aTwM[ii].m_fX;
		m_aTwM[ii].m_fYo = m_aTwM[ii].m_fY;
		m_aTwM[ii].m_fZo = m_aTwM[ii].m_fZ;

		m_aBwM[ii].m_fVx = 0.0;
		m_aBwM[ii].m_fVy = 0.0;
		m_aBwM[ii].m_fVz = 0.0;
	
		m_aTwM[ii].m_fVx = 0.0;
		m_aTwM[ii].m_fVy = 0.0;
		m_aTwM[ii].m_fVz = 0.0;


	
	}



	SortGridInit();
//	SortWallMolecules();
	//OnRun();


}
