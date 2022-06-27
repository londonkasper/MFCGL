// PropDlg.cpp : ���� �����Դϴ�.
//

#include "stdafx.h"
#include "MFCGL.h"
#include "PropDlg.h"


// CPropDlg ��ȭ �����Դϴ�.

IMPLEMENT_DYNAMIC(CPropDlg, CDialog)

CPropDlg::CPropDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CPropDlg::IDD, pParent)
{

}

CPropDlg::~CPropDlg()
{
}

void CPropDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CPropDlg, CDialog)
	ON_BN_CLICKED(IDOK, &CPropDlg::OnBnClickedOk)
	ON_EN_CHANGE(IDC_m_nroop, &CPropDlg::OnEnChangemnroop)
	ON_EN_CHANGE(IDC_m_gas, &CPropDlg::OnEnChangemgas)
END_MESSAGE_MAP()


// CPropDlg �޽��� ó�����Դϴ�.

void CPropDlg::OnBnClickedOk()
{


	m_nwMs = GetDlgItemInt(IDC_m_nwMs);
	m_nMs = GetDlgItemInt(IDC_m_niMs);
	m_fMrad = GetDlgItemInt(IDC_m_fMrad);
	m_fMdist = GetDlgItemInt(IDC_m_fMdist);
	//m_fH = GetDlgItemInt(IDC_m_fH);
	m_fH2 = GetDlgItemInt(IDC_m_fH2);
	m_fW2 = GetDlgItemInt(IDC_m_fW2);
	m_fZ2 = GetDlgItemInt(IDC_m_fZ2);
	m_nWtype = GetDlgItemInt(IDC_m_Wtype);
	m_ncutoff = GetDlgItemInt(IDC_m_cutoff);
	m_nwallinteraction = GetDlgItemInt(IDC_m_wallinteraction);
	m_ndamping = GetDlgItemInt(IDC_m_damping);
	m_nroop = GetDlgItemInt(IDC_m_nroop);
	m_nmilisec = GetDlgItemInt(IDC_m_nmilisec);
	m_dt = GetDlgItemInt(IDC_m_dt);
	m_targett = GetDlgItemInt(IDC_m_targett);
	m_tt = GetDlgItemInt(IDC_m_tt);
	m_bt = GetDlgItemInt(IDC_m_bt);
	m_spring = GetDlgItemInt(IDC_m_spring);

	char   *stopstring;
	CString str1;
	GetDlgItemText(IDC_m_istrt,str1);
    m_istrt= strtod(str1,&stopstring);

	
	CString str3;
	GetDlgItemText(IDC_m_gas,str3);
    m_acc= strtod(str3,&stopstring);


	CString str5;
	GetDlgItemText(IDC_m_max,str5);
    m_max= strtod(str5,&stopstring);

	OnOK();
}

BOOL CPropDlg::OnInitDialog() 
{
   CDialog::OnInitDialog();


	SetDlgItemInt(IDC_m_nwMs,m_nwMs);
	SetDlgItemInt(IDC_m_niMs,m_nMs);
	SetDlgItemInt(IDC_m_fMrad,m_fMrad);
	SetDlgItemInt(IDC_m_fMdist,m_fMdist);
	//SetDlgItemInt(IDC_m_fH,m_fH);
	SetDlgItemInt(IDC_m_fH2,m_fH2);
	SetDlgItemInt(IDC_m_fW2,m_fW2);
	SetDlgItemInt(IDC_m_fZ2,m_fZ2);
	SetDlgItemInt(IDC_m_Wtype,m_nWtype);
	SetDlgItemInt(IDC_m_cutoff,m_ncutoff);
	SetDlgItemInt(IDC_m_wallinteraction,m_nwallinteraction);
	SetDlgItemInt(IDC_m_damping,m_ndamping);
	SetDlgItemInt(IDC_m_nroop,m_nroop);
	SetDlgItemInt(IDC_m_nmilisec,m_nmilisec);
	SetDlgItemInt(IDC_m_dt,m_dt);
	SetDlgItemInt(IDC_m_targett,m_targett);
	SetDlgItemInt(IDC_m_tt,m_tt);
	SetDlgItemInt(IDC_m_bt,m_bt);
	SetDlgItemInt(IDC_m_spring,m_spring);


	CString str2;
	str2.Format("%lf",m_istrt);
	SetDlgItemTextA(IDC_m_istrt,str2);

	CString str4;
	str4.Format("%lf",m_acc);
	SetDlgItemTextA(IDC_m_acc,str4);

	CString str6;
	str6.Format("%lf",m_max);
	SetDlgItemTextA(IDC_m_max,str6);

   return TRUE; 
}

void CPropDlg::OnEnChangemnroop()
{
	// TODO:  RICHEDIT ��Ʈ���� ���, �� ��Ʈ����
	// CDialog::OnInitDialog() �Լ��� �������ϰ�  ����ũ�� OR �����Ͽ� ������
	// ENM_CHANGE �÷��׸� �����Ͽ� CRichEditCtrl().SetEventMask()�� ȣ���ؾ߸�
	// �ش� �˸� �޽����� �����ϴ�.

	// TODO:  ���⿡ ��Ʈ�� �˸� ó���� �ڵ带 �߰��մϴ�.
}

void CPropDlg::OnEnChangemgas()
{
	// TODO:  RICHEDIT ��Ʈ���� ���, �� ��Ʈ����
	// CDialog::OnInitDialog() �Լ��� �������ϰ�  ����ũ�� OR �����Ͽ� ������
	// ENM_CHANGE �÷��׸� �����Ͽ� CRichEditCtrl().SetEventMask()�� ȣ���ؾ߸�
	// �ش� �˸� �޽����� �����ϴ�.

	// TODO:  ���⿡ ��Ʈ�� �˸� ó���� �ڵ带 �߰��մϴ�.
}
