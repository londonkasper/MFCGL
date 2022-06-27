#pragma once


// CSDlg 대화 상자입니다.

class CSDlg : public CDialog
{
	DECLARE_DYNAMIC(CSDlg)

public:
	CSDlg(CWnd* pParent = NULL);   // 표준 생성자입니다.
	virtual ~CSDlg();

	int m_nnz;
	int m_nnx;


// 대화 상자 데이터입니다.
	enum { IDD = IDD_DIALOG3 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 지원입니다.

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedOk();
};
