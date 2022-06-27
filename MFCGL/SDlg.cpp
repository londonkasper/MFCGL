// SDlg.cpp : 구현 파일입니다.
//

#include "stdafx.h"
#include "MFCGL.h"
#include "SDlg.h"



// CSDlg 대화 상자입니다.

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


// CSDlg 메시지 처리기입니다.

BOOL CSDlg::OnInitDialog()
{


	CDialog::OnInitDialog();

	SetDlgItemInt(IDC_EDIT1,m_nnx);
	SetDlgItemInt(IDC_EDIT3,m_nnz);


	return TRUE;  // return TRUE unless you set the focus to a control
	// 예외: OCX 속성 페이지는 FALSE를 반환해야 합니다.
}

void CSDlg::OnBnClickedOk()
{
	m_nnx = GetDlgItemInt(IDC_EDIT1);
	m_nnz = GetDlgItemInt(IDC_EDIT3);

	OnOK();
}
