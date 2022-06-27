// LatticeDlg.cpp : 구현 파일입니다.
//

#include "stdafx.h"
#include "MFCGL.h"
#include "LatticeDlg.h"


// CLatticeDlg 대화 상자입니다.

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


// CLatticeDlg 메시지 처리기입니다.

void CLatticeDlg::OnEnChangeEdit1()
{
	// TODO:  RICHEDIT 컨트롤인 경우, 이 컨트롤은
	// CDialog::OnInitDialog() 함수를 재지정하고  마스크에 OR 연산하여 설정된
	// ENM_CHANGE 플래그를 지정하여 CRichEditCtrl().SetEventMask()를 호출해야만
	// 해당 알림 메시지를 보냅니다.

	// TODO:  여기에 컨트롤 알림 처리기 코드를 추가합니다.
}

void CLatticeDlg::OnBnClickedOk()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.

	
	m_nx = GetDlgItemInt(IDC_EDIT1);
	m_ny = GetDlgItemInt(IDC_EDIT2);
	m_nz = GetDlgItemInt(IDC_EDIT3);

	OnOK();
}
