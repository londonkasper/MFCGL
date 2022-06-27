#include "StdAfx.h"
#include "Molecular.h"
#include "MFCGL.h"

CMolecular::CMolecular(void)
{

	m_fX = 0.0;
	m_fY = 0.0;
	m_fZ = 0.0;

	m_fXo = 0.0;
	m_fYo = 0.0;
	m_fZo = 0.0;

	m_fVx = 0.0;
	m_fVy = 0.0;
	m_fVz = 0.0;

	m_fax = 0.0;
	m_fay = 0.0;
	m_faz = 0.0;

	m_fMmass = 1.0;

	m_fPot = 0.0;
	m_fForce = 0.0;
	m_fStress = 0.0;
	m_fVirialxy = 0.0;
	m_dh = 0.0;

}

CMolecular::~CMolecular(void)
{
}
