// LatticeDlg.cpp : ���� �����Դϴ�.
//

#include "stdafx.h"
#include "MFCGL.h"
#include "LatticeDlg.h"


// CLatticeDlg ��ȭ �����Դϴ�.

IMPLEMENT_DYNAMIC(CLatticeDlg, CDialog)

CLatticeDlg::CLatticeDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CLatticeDlg::IDD, pParent)
{

}

BOOL CLatticeDlg::OnInitDialog() 
{
   CDialog::OnInitDialog();


	SetDlgItemInt(IDC_EDIT1,m_nx);
	SetDlgItemInt(IDC_EDIT2,m_ny);
	SetDlgItemInt(IDC_EDIT3,m_nz);

	return TRUE; 

}

CLatticeDlg::~CLatticeDlg()
{
}

void CLatticeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CLatticeDlg, CDialog)
	ON_EN_CHANGE(IDC_EDIT1, &CLatticeDlg::OnEnChangeEdit1)
	ON_BN_CLICKED(IDOK, &CLatticeDlg::OnBnClickedOk)
END_MESSAGE_MAP()


// CLatticeDlg �޽��� ó�����Դϴ�.

void CLatticeDlg::OnEnChangeEdit1()
{
	// TODO:  RICHEDIT ��Ʈ���� ���, �� ��Ʈ����
	// CDialog::OnInitDialog() �Լ��� �������ϰ�  ����ũ�� OR �����Ͽ� ������
	// ENM_CHANGE �÷��׸� �����Ͽ� CRichEditCtrl().SetEventMask()�� ȣ���ؾ߸�
	// �ش� �˸� �޽����� �����ϴ�.

	// TODO:  ���⿡ ��Ʈ�� �˸� ó���� �ڵ带 �߰��մϴ�.
}

void CLatticeDlg::OnBnClickedOk()
{
	// TODO: ���⿡ ��Ʈ�� �˸� ó���� �ڵ带 �߰��մϴ�.

	
	m_nx = GetDlgItemInt(IDC_EDIT1);
	m_ny = GetDlgItemInt(IDC_EDIT2);
	m_nz = GetDlgItemInt(IDC_EDIT3);

	OnOK();
}
