// SDlg.cpp : ���� �����Դϴ�.
//

#include "stdafx.h"
#include "MFCGL.h"
#include "SDlg.h"



// CSDlg ��ȭ �����Դϴ�.

IMPLEMENT_DYNAMIC(CSDlg, CDialog)

CSDlg::CSDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CSDlg::IDD, pParent)
{

}

CSDlg::~CSDlg()
{
}

void CSDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CSDlg, CDialog)
	ON_BN_CLICKED(IDOK, &CSDlg::OnBnClickedOk)
END_MESSAGE_MAP()


// CSDlg �޽��� ó�����Դϴ�.

BOOL CSDlg::OnInitDialog()
{


	CDialog::OnInitDialog();

	SetDlgItemInt(IDC_EDIT1,m_nnx);
	SetDlgItemInt(IDC_EDIT3,m_nnz);


	return TRUE;  // return TRUE unless you set the focus to a control
	// ����: OCX �Ӽ� �������� FALSE�� ��ȯ�ؾ� �մϴ�.
}

void CSDlg::OnBnClickedOk()
{
	m_nnx = GetDlgItemInt(IDC_EDIT1);
	m_nnz = GetDlgItemInt(IDC_EDIT3);

	OnOK();
}
