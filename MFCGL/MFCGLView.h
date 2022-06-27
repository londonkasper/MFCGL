// MFCGLView.h : interface of the CMFCGLView class
//


#include "Molecular.h"
#include "PropDlg.h"
#include "LatticeDlg.h"
#include "SDlg.h"
//#include <fstream.h>

#pragma once


class CMFCGLView : public CView
{
protected: // create from serialization only
	CMFCGLView();
	DECLARE_DYNCREATE(CMFCGLView)

// Attributes
public:
	CMFCGLDoc* GetDocument() const;



	CStdioFile m_file;  //Export file
	CStdioFile m_rfile; //Read file

	CStdioFile m_ygridfile; //Read file  //data collect from ygrid
	CStdioFile m_binfile; //Read file	 //data collect from bins
	CStdioFile m_heatfluxfile;


	CMolecular *m_aTwM;
	CMolecular *m_aBwM;

	CMolecular *m_aTwM_OMP;
	CMolecular *m_aBwM_OMP;

	CMolecular *m_aTwM_old;
	CMolecular *m_aBwM_old;
	CMolecular *m_aTwM_orig;
	CMolecular *m_aBwM_orig;

	CMolecular *m_aM;
	CMolecular *m_aM_OMP;

	CMolecular *m_aMVel;

////MB////////
	int visual, Nfolder;
	CString folder;
	char bufferfolder[20];
////MB///////

	struct InterVector            // interaction force vector
	{
		double fsignx;
		double fsigny;
		double fsignz;
		double fdist,fang,fForce,fPot;
		double accx,accy,accz;
	};      

	InterVector m_vec;

	struct Cell_Struct            // Declare PERSON struct type
	{
	int      NM;
	int      N[64]; 
	double   x[64];              // Declare member types
	double   y[64];              // Declare member types
	double   z[64];              // Declare member types
	double   vx[64];              // Declare member types
	double   vy[64];              // Declare member types
	double   vz[64];              // Declare member types
	int mtype[64];                   // 1: fluid, 2: wall, 3: image
	};      


	Cell_Struct Cell[29][27][11];
/*
	Cell_Struct Bin[50][60];

	Cell_Struct ImagR[1500];
	Cell_Struct ImagL[1500];

	Cell_Struct WallT[1500];
	Cell_Struct WallB[1500];

	Cell_Struct WimgTR[1];
	Cell_Struct WimgTL[1];
	Cell_Struct WimgBR[1];
	Cell_Struct WimgBL[1];
*/	
	int m_aNMVel[200];  //For Number of particle represneting vel or density
	
	CString m_strTemp;
	CString m_strP;
	CString m_strNM;
	CString m_strTWallTemp;
	CString m_strBWallTemp;


	CString strG[10];
	CString strGDen[10];
	double dG[10];   //10 section temperature
	double denCount[10]; ////10 section number density

	
	//Coordinate variable
	double m_xRotate;
	double m_yRotate;
	double m_zRotate;
	double m_xTranslate;
	double m_yTranslate;
	double m_zTranslate;
	double m_dautocolor;
	double view_angle;

	//Input dialog
	CPropDlg m_Propdlg; 
	CLatticeDlg m_LatticeDlg; 
	CSDlg m_SDlg; 

	//Propertys of the simulation
	int m_nwMs; //Number of Molecules
	int m_nMs;  //Number of Molecules
	double m_fMrad;  //Radius of the sphere
	double m_fMdist; //Distance between wall molecule
	double m_fdt; //Time step
	double m_dcutdist;
	double m_dEwall,m_dEgas;
	double m_dMass_wall,m_dMass_gas;
	double m_AvogN;
	double m_UGasN;
	double m_BoltzN;
	double m_dGravP;
	double m_dSigma;
	double m_cutpot;
	double m_fH;  //Wall Separation
	double m_fW;  //Channel Length
	double m_fZ;  //Depth(Z)
	double m_dampcoef;
	double m_Ttarget; // Target Temperature
	double m_cellsize;
	double m_ygsize;
	double m_bxsize;
	double m_bysize;
	double d_vibw; //Wall vibration frequency
	double d_numden;


	int m_fH2;  
	int m_fW2;  
	int m_fZ2;  

	int m_nz;
	int m_nx;
	int m_ny;


	int m_nnz;
	int m_nnx;
	int m_nny;

	double m_istrt;
	double m_acc;
	double m_max;
	double tau;

	double d_Twalltemp; //Wall vibration frequency
	double d_Bwalltemp; //Wall vibration frequency

	double d_Twalltemp1; //Wall vibration frequency
	double d_Bwalltemp1; //Wall vibration frequency
	double d_Twalltemp2; //Wall vibration frequency
	double d_Bwalltemp2; //Wall vibration frequency
	double d_Twalltemp3; //Wall vibration frequency
	double d_Bwalltemp3; //Wall vibration frequency

	double wallTempT,wallTempB; //target wall temp

	double d_K_M; // K-divided by m for stiffness

	double m_massofargon;

	double d_ntimestep;
	double d_totalLoop;

	//for unit conversion
	double joule_to_gnmps;
	double m_to_nm,s_to_ps,kg_to_g;

	double m_heatfluxx,m_heatfluxy,m_heatfluxz,m_heatsumx,m_heatsumy,m_heatsumz;
	double d_Jx,d_Jy,d_Jz;  //Heat flux
	double d_J[100000][3];
	int  	n_Jcount,n_Jtotal;

	double m_V;

	double m_V1[1000]; //for wall vibration
	double m_V2[1000];

	//For Nose-Hoover Thermostat;
	double m_zeta;
	double m_zeta_dot;
	double m_tau;
	double m_Q;

	double m_Umax;  //velocity of top and bottom wall
	double m_gamma; //shear rate of the fluid
	double m_gamma2; //for the pressure driven flow


	double m_T;
	double m_P;
	double m_Kin;
	double m_Pot;
	double m_Potx,m_Poty;
	double m_Tot;

	double m_Tsc;
	double m_Kinsc;
	double m_Potsc;
	double m_Totsc;

	double m_Tsca[400];
	double m_Kinsca[400];
	double m_Potsca[400];
	double m_Totsca[400];

		double m_Potsca1[400];
		double m_Potsca2[400];
		double m_Potsca3[400];
		double m_Potsca4[400];

		double m_Kinsca1[400];
		double m_Kinsca2[400];
		double m_Kinsca3[400];
		double m_Kinsca4[400];

	int m_Gcount;

	double datacollect[2000];
	double datacollect1[2000];
	double datacollect2[2000];
	double bin_datacollect[2000];


	double m_Kb; 
	double m_rf;   //Reverse plot


	//Flag for the mode [0:Off 1:On]
	int m_nMode;
	int m_nShow_molecule;
	int m_nWire;
	int m_nPeriodic;
	int m_nSp;
	int m_nWtype;
	int m_n3d;
	int m_ncutoff;
	int m_nwallinteraction;
	int m_ndamping;
	int m_time;
	int m_nNrow;
	int m_nWrange;
	int m_nroop;
	int m_nmilisec;
	int m_nWallDown;
	int m_spring;

	//For data collection and time averaging!
	double m_data_count5[500];
	double m_data_count10[500];
	double m_data_count50[500];
	double m_data_count100[500];
	double m_data_count500[500];

	double m_data_temp5[500];
	double m_data_temp10[500];
	double m_data_temp50[500];
	double m_data_temp100[500];
	double m_data_temp500[500];

	double m_data_vel5[500];
	double m_data_vel10[500];
	double m_data_vel50[500];
	double m_data_vel100[500];
	double m_data_vel500[500];

	double m_data_collectT5[500];
	double m_data_collectT10[500];
	double m_data_collectT50[500];
	double m_data_collectT100[500];
	double m_data_collectT500[500];

	double m_data_collectD5[500];
	double m_data_collectD10[500];
	double m_data_collectD50[500];
	double m_data_collectD100[500];
	double m_data_collectD500[500];

	double m_data_collectV5[500];
	double m_data_collectV10[500];
	double m_data_collectV50[500];
	double m_data_collectV100[500];
	double m_data_collectV500[500];


	double m_data_stress5[500];
	double m_data_stress10[500];
	double m_data_stress100[500];

	double m_data_virial5[500];
	double m_data_virial10[500];
	double m_data_virial100[500];

	double m_data_collectstress5[500];
	double m_data_collectstress10[500];
	double m_data_collectstress100[500];

	double m_data_collectvirial5[500];
	double m_data_collectvirial10[500];
	double m_data_collectvirial100[500];


	int m_startstep1;
	int m_endstep1;
	int m_averagingstep1;

	int m_startstep2;
	int m_endstep2;
	int m_averagingstep2;
	
	int m_startstep3;
	int m_endstep3;
	int m_averagingstep3;
	
	int m_startstep4;
	int m_endstep4;
	int m_averagingstep4;
	
	int m_startstep5;
	int m_endstep5;
	int m_averagingstep5;
	//For data collection and time averaging!
	

	CString m_str_NumberofStep,m_str_totalLoop,m_str_heatfluxx,m_str_heatfluxy,m_str_heatfluxz;
	CString BT1,BT2,BT3,TT1,TT2,TT3;

/////////////////////////////


	int nRotate;
	int iCPU;
	double glcrossx;
	double glcrossy;

	//for calctime measure
	clock_t start, finish;
	double  duration;

	
	
	//clock_t start, finish;
	//double  duration;
	//start = clock();
	//finish = clock();
	//duration = (double)(finish - start) / CLOCKS_PER_SEC;






// Operations
public:



// Overrides
public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);



// Implementation
public:

	void TmpFileSave();
	void AssignImage(int ai,int bi,int ci,int ao,int bo,int co, double xx, double yy, double zz);
	void AssignPerBC();
	void SortWallMolecules();
	void SortGridInit();
	void CollectData();
	void CollectEnergy();

	void DriveFlow();
	void GetCellInteraction();
	void GetReactionForce(double xi,double xj, double yi,double yj, double zi,double zj,double e);
	void CollectAndSaveFile();
	void SetInitialCondition();

////MB////////
	void visualclose();
////MB///////

	virtual ~CMFCGLView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

	CDC* m_pDC;   // 디바이스 컨텍스트(DC)
	HGLRC m_hRC;  // 렌더링 컨텍스트(RC)
//	HDC m_hDC;
	void DrawScene(void);
	BOOL SetupPixelFormat(PIXELFORMATDESCRIPTOR* pPFD=0);
	GLuint base;
	GLvoid BuildFont(GLvoid);
	int InitGL(GLvoid);
	GLvoid glPrint(const char *text);
	GLvoid KillFont(GLvoid);

	void GLResize(int cx, int cy);
	void OMPTest(void);

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:


	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);


	afx_msg void OnMove();
	afx_msg void OnShowmolecule();
	afx_msg void OnNewsim();
	afx_msg void OnWire();
	afx_msg void OnRun();
	afx_msg void OnShear();
	afx_msg void OnPeriodic();
	afx_msg void OnPressure();
	afx_msg void OnInit();
	afx_msg void OnPropdlg();

	afx_msg void OnTimer(UINT_PTR nIDEvent);

	afx_msg void OnAutoTimer();
	afx_msg void OnDamp1();
	afx_msg void OnDampdn();
	afx_msg void OnDampup();

	afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);




	afx_msg void OnAutoRun();
	afx_msg void OnAutoStop();
	afx_msg void OnAutoSinglecpu();
	afx_msg void OnAutoAddcpu();
	afx_msg void OnAutoMaxcpu();
	afx_msg void OnAutoRemovecpu();



	afx_msg void OnFileOpen();
	afx_msg void OnFileSave();
	afx_msg void OnWalldown();
	afx_msg void OnReverseflot();
	afx_msg void OnBoltzmannVel();
	afx_msg void OnOpenwall();
	afx_msg void OnMsave();
	afx_msg void OnPost();
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnLattice();
	afx_msg void OnS();
};

#ifndef _DEBUG  // debug version in MFCGLView.cpp
inline CMFCGLDoc* CMFCGLView::GetDocument() const
   { return reinterpret_cast<CMFCGLDoc*>(m_pDocument); }
#endif

