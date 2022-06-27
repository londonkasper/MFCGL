// MFCGLDoc.cpp : implementation of the CMFCGLDoc class
//

#include "stdafx.h"
#include "MFCGL.h"

#include "MFCGLDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMFCGLDoc

IMPLEMENT_DYNCREATE(CMFCGLDoc, CDocument)

BEGIN_MESSAGE_MAP(CMFCGLDoc, CDocument)
END_MESSAGE_MAP()


// CMFCGLDoc construction/destruction

CMFCGLDoc::CMFCGLDoc()
{
	// TODO: add one-time construction code here

}

CMFCGLDoc::~CMFCGLDoc()
{
}

BOOL CMFCGLDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}




// CMFCGLDoc serialization

void CMFCGLDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}


// CMFCGLDoc diagnostics

#ifdef _DEBUG
void CMFCGLDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CMFCGLDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CMFCGLDoc commands
