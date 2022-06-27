// MFCGL.h : main header file for the MFCGL application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols


// CMFCGLApp:
// See MFCGL.cpp for the implementation of this class
//

class CMFCGLApp : public CWinApp
{
public:
	CMFCGLApp();


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CMFCGLApp theApp;