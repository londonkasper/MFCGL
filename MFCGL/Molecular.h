#pragma once

class CMolecular
{
public:
	CMolecular(void);
	~CMolecular(void);


	double m_fX,m_fY,m_fZ;    //Position
	double m_fXo,m_fYo,m_fZo; //Original Position

	double m_fVx,m_fVy,m_fVz; //Velocity

	double m_fax,m_fay,m_faz;

	double m_fMmass;

	double m_fPot;
	double m_fForce;
	double m_fStress;
	double m_fVirialxy;
	double m_dh;



};
