#pragma once


// CLatticeDlg ��ȭ �����Դϴ�.

class CLatticeDlg : public CDialog
{
	DECLARE_DYNAMIC(CLatticeDlg)

	virtual BOOL OnInitDialog();

public:
	CLatticeDlg(CWnd* pParent = NULL);   // ǥ�� �������Դϴ�.
	virtual ~CLatticeDlg();

	int m_nz;
	int m_nx;
	int m_ny;


// ��ȭ ���� �������Դϴ�.
	enum { IDD = IDD_DIALOG2 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV �����Դϴ�.

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedOk();
};
