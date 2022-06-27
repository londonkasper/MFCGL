#pragma once


// CLatticeDlg 대화 상자입니다.

class CLatticeDlg : public CDialog
{
	DECLARE_DYNAMIC(CLatticeDlg)

	virtual BOOL OnInitDialog();

public:
	CLatticeDlg(CWnd* pParent = NULL);   // 표준 생성자입니다.
	virtual ~CLatticeDlg();

	int m_nz;
	int m_nx;
	int m_ny;


// 대화 상자 데이터입니다.
	enum { IDD = IDD_DIALOG2 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 지원입니다.

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedOk();
};
