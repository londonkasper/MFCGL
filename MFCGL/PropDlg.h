#pragma once


// CPropDlg ��ȭ �����Դϴ�.

class CPropDlg : public CDialog
{
	DECLARE_DYNAMIC(CPropDlg)

public:
	CPropDlg(CWnd* pParent = NULL);   // ǥ�� �������Դϴ�.
	virtual ~CPropDlg();

	virtual BOOL OnInitDialog();


	int m_nwMs;           //Number of Wall Molecules
	int m_nMs;         //Number of Inner Molecules

	int m_fMrad;         //Radius of Sphere
	int m_fMdist;        //Distance between Molecules 
	int m_fH;   //Width of the Channel
	int m_fW;
	int m_fZ;
	int m_nWtype;            //Wall type 0: BCC  1:FCC
	int m_ncutoff;
	int m_nwallinteraction;
	int m_ndamping;
	int m_nroop;
	int m_nmilisec;
	int m_dt;
	int m_targett;
	int m_tt;
	int m_bt;
	int m_spring;
	int m_fH2;  
	int m_fW2;  
	int m_fZ2;  
	double m_istrt;
	double m_acc;
	double m_max;




// ��ȭ ���� �������Դϴ�.
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV �����Դϴ�.

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnEnChangemnroop();
	afx_msg void OnEnChangemgas();
};
