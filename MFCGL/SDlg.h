#pragma once


// CSDlg ��ȭ �����Դϴ�.

class CSDlg : public CDialog
{
	DECLARE_DYNAMIC(CSDlg)

public:
	CSDlg(CWnd* pParent = NULL);   // ǥ�� �������Դϴ�.
	virtual ~CSDlg();

	int m_nnz;
	int m_nnx;


// ��ȭ ���� �������Դϴ�.
	enum { IDD = IDD_DIALOG3 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV �����Դϴ�.

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedOk();
};
