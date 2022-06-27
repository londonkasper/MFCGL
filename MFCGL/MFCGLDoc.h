// MFCGLDoc.h : interface of the CMFCGLDoc class
//


#pragma once


class CMFCGLDoc : public CDocument
{
protected: // create from serialization only
	CMFCGLDoc();
	DECLARE_DYNCREATE(CMFCGLDoc)

// Attributes
public:

// Operations
public:

// Overrides
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CMFCGLDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
};


