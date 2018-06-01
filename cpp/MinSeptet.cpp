#include<vector>
#include<cstring>
#include<algorithm>
#include <iostream>

#include "CSplitSys.h"
#include "auxilary.h" 

using namespace std;

// calculate w(321|7654):72 possible formulas 

 double CSplitSys::MinSeptet(unsigned int gLeft[3], unsigned int gRight[4])
{  
    CSextet sexTemp;
    unsigned int a, b, c, d, e, f, g;
    a = gLeft[2]; //1
    b = gLeft[1]; //2
    c = gLeft[0]; //3
    d = gRight[3]; //4
    e = gRight[2]; //5
    f = gRight[1]; //6
    g = gRight[0]; //7
    //all possible left block or right block may be used in the formulation
    unsigned int B321[3] = {c, b, a};
    unsigned int B421[3] = {d, b, a};
    unsigned int B431[3] = {d, c, a};
    unsigned int B521[3] = {e, b, a};
    unsigned int B531[3] = {e, c, a};
    unsigned int B541[3] = {e, d, a};
    unsigned int B621[3] = {f, b, a};
    unsigned int B631[3] = {f, c, a};
    unsigned int B641[3] = {f, d, a};
    unsigned int B651[3] = {f, e, a};
    unsigned int B721[3] = {g, b, a};
    unsigned int B731[3] = {g, c, a};
    unsigned int B741[3] = {g, d, a};
    unsigned int B751[3] = {g, e, a};
    unsigned int B761[3] = {g, f, a};
    
    unsigned int B432[3] = {d, c, b};
    unsigned int B532[3] = {e, c, b};
    unsigned int B542[3] = {e, d, b};
    unsigned int B632[3] = {f, c, b};
    unsigned int B642[3] = {f, d, b};
    unsigned int B652[3] = {f, e, b};
    unsigned int B732[3] = {g, c, b};
    unsigned int B742[3] = {g, d, b};
    unsigned int B752[3] = {g, e, b};
    unsigned int B762[3] = {g, f, b};
    
    unsigned int B543[3] = {e, d, c};
    unsigned int B643[3] = {f, d, c};
    unsigned int B653[3] = {f, e, c};
    unsigned int B743[3] = {g, d, c};
    unsigned int B753[3] = {g, e, c};
    unsigned int B763[3] = {g, f, c};
    
    unsigned int B654[3] = {f, e, d};
    unsigned int B754[3] = {g, e, d};
    unsigned int B764[3] = {g, f, d};
    
    unsigned int B765[3] = {g, f, e};
    
    //calculate all the possible weight may be used in the formulas
    sexTemp.SetLabel(B654, B321);
    double w654321 = m_gSextet[sexTemp.GetIndex()].GetWeight();

    sexTemp.SetLabel(B721, B654);
    double w721654 = m_gSextet[sexTemp.GetIndex()].GetWeight();
   
    sexTemp.SetLabel(B721, B543);
    double w721543 = m_gSextet[sexTemp.GetIndex()].GetWeight();
   
    sexTemp.SetLabel(B761, B543);
    double w761543 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B761, B432);
    double w761432 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B765, B432);
    double w765432 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B765, B321);
    double w765321 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B761, B532);
    double w761532 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B764, B532);
    double w764532 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B764, B321);
    double w764321 = m_gSextet[sexTemp.GetIndex()].GetWeight();
   
    sexTemp.SetLabel(B762, B421);
    double w762421 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B765, B431);
    double w765431 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B762, B531);
    double w762531 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B762, B543);
    double w762543 = m_gSextet[sexTemp.GetIndex()].GetWeight();
   
    sexTemp.SetLabel(B764, B531);
    double w764531 = m_gSextet[sexTemp.GetIndex()].GetWeight();
   
    sexTemp.SetLabel(B721, B643);
    double w721643 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B751, B643);
    double w751643 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B751, B432);
    double w751432 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B751, B632);
    double w751632 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B754, B632);
    double w754632 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B754, B321);
    double w754321 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B752, B643);
    double w752643 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B752, B431);
    double w752431 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B752, B631);
    double w752631 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B754, B631);
    double w754631 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B721, B653);
    double w721653 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B741, B653);
    double w741653 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B741, B532);
    double w741532 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B741, B632);
    double w741632 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B742, B653);
    double w742653 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B742, B531);
    double w742531 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B742, B631);
    double w742631 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B731, B654);
    double w731654 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B731, B542);
    double w731542 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B761, B542);
    double w761542 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B763, B542);
    double w763542 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B763, B421);
    double w763421 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B765, B421);
    double w765421 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B763, B521);
    double w763521 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B764, B521);
    double w764521 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B751, B642);
    double w751642 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B753, B642);
    double w753642 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B753, B421);
    double w753421 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B753, B621);
    double w753621 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B754, B621);
    double w754621 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B731, B652);
    double w731652 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B741, B652);
    double w741652 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B743, B652);
    double w743652 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B743, B521);
    double w743521 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B743, B621);
    double w743621 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B732, B654);
    double w732654 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B732, B541);
    double w732541 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B762, B541);
    double w762541 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B763, B541);
    double w763541 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B732, B641);
    double w732641 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B752, B641);
    double w752641 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B753, B641);
    double w753641 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B732, B651);
    double w732651 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B742, B651);
    double w742651 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B743, B651);
    double w743651 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B621, B543);
    double w621543 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B651, B432);
    double w651432 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B652, B431);
    double w652431 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B641, B532);
    double w641532 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B642, B531);
    double w642531 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B631, B542);
    double w631542 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B653, B421);
    double w653421 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B643, B521);
    double w643521 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B632, B541);
    double w632541 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B731, B642);
    double w731642 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    sexTemp.SetLabel(B762, B431);
    double w762431 = m_gSextet[sexTemp.GetIndex()].GetWeight();
    
    //calculate the formulations
    //1:w(123|456)-w(127|456)+w(127|345)-w(167|345)+w(167|234)-w(567|234)+w(567|123)
    double dTemp = 0;
    dTemp = w654321 - w721654 + w721543 - w761543 + w761432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    
    double dWei =dTemp;
    
    //2:w(123|456)-w(127|456)+w(127|345)-w(167|345)+w(167|235)-w(467|235)+w(467|123)
    dTemp = w654321 - w721654 + w721543 - w761543 + w761532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //3:w(123|456)-w(127|456)+w(127|345)-w(267|345)+w(267|134)-w(567|134)+w(567|123)
    dTemp = w654321 - w721654 + w721543 - w762543 + w762431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //4:w(123|456)-w(127|456)+w(127|345)-w(267|345)+w(267|135)-w(467|135)+w(467|123)
    dTemp = w654321 - w721654 + w721543 - w762543 + w762531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //5:w(123|456)-w(127|456)+w(127|346)-w(157|346)+w(157|234)-w(567|234)+w(567|123)
    dTemp = w654321 - w721654 + w721643 - w751643 + w751432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //6:w(123|456)-w(127|456)+w(127|346)-w(157|346)+w(157|236)-w(457|236)+w(457|123)
    dTemp = w654321 - w721654 + w721643 - w751643 + w751632 - w754632 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //7:w(123|456)-w(127|456)+w(127|346)-w(257|346)+w(257|134)-w(567|134)+w(567|123)
    dTemp = w654321 - w721654 + w721643 - w752643 + w752431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //8:w(123|456)-w(127|456)+w(127|346)-w(257|346)+w(257|136)-w(457|136)+w(457|123)
    dTemp = w654321 - w721654 + w721643 - w752643 + w752631 - w754631 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //9:w(123|456)-w(127|456)+w(127|356)-w(147|356)+w(147|235)-w(467|235)+w(467|123)
    dTemp = w654321 - w721654 + w721653 - w741653 + w741532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //10:w(123|456)-w(127|456)+w(127|356)-w(147|356)+w(147|236)-w(457|236)+w(457|123)
    dTemp = w654321 - w721654 + w721653 - w741653 + w741632 - w754632 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //11:w(123|456)-w(127|456)+w(127|356)-w(247|356)+w(247|135)-w(467|135)+w(467|123)
    dTemp = w654321 - w721654 + w721653 - w742653 + w742531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //12:w(123|456)-w(127|456)+w(127|356)-w(247|356)+w(247|136)-w(457|136)+w(457|123)
    dTemp = w654321 - w721654 + w721653 - w742653 + w742631 - w754631 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //13:w(123|456)-w(137|456)+w(137|245)-w(167|245)+w(167|234)-w(567|234)+w(567|123)
    dTemp = w654321 - w731654 + w731542 - w761542 + w761432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //14:w(123|456)-w(137|456)+w(137|245)-w(167|245)+w(167|235)-w(467|235)+w(467|123)
    dTemp = w654321 - w731654 + w731542 - w761542 + w761532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //15:w(123|456)-w(137|456)+w(137|245)-w(367|245)+w(367|124)-w(567|124)+w(567|123)
    dTemp = w654321 - w731654 + w731542 - w763542 + w763421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //16:w(123|456)-w(137|456)+w(137|245)-w(367|245)+w(367|125)-w(467|125)+w(467|123)
    dTemp = w654321 - w731654 + w731542 - w763542 + w763521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //17:w(123|456)-w(137|456)+w(137|246)-w(157|246)+w(157|234)-w(567|234)+w(567|123)
    dTemp = w654321 - w731654 + w731642 - w751642 + w751432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //18:w(123|456)-w(137|456)+w(137|246)-w(157|246)+w(157|236)-w(457|236)+w(457|123)
    dTemp = w654321 - w731654 + w731642 - w751642 + w751632 - w754632 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //19:w(123|456)-w(137|456)+w(137|246)-w(357|246)+w(357|124)-w(567|124)+w(567|123)
    dTemp = w654321 - w731654 + w731642 - w753642 + w753421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //20:w(123|456)-w(137|456)+w(137|246)-w(357|246)+w(357|126)-w(457|126)+w(457|123)
    dTemp = w654321 - w731654 + w731642 - w753642 + w753621 - w754621 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //21:w(123|456)-w(137|456)+w(137|256)-w(147|256)+w(147|235)-w(467|235)+w(467|123)
    dTemp = w654321 - w731654 + w731652 - w741652 + w741532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //22:w(123|456)-w(137|456)+w(137|256)-w(147|256)+w(147|236)-w(457|236)+w(457|123)
    dTemp = w654321 - w731654 + w731652 - w741652 + w741632 - w754632 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //23:w(123|456)-w(137|456)+w(137|256)-w(347|256)+w(347|125)-w(467|125)+w(467|123)
    dTemp = w654321 - w731654 + w731652 - w743652 + w743521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //24:w(123|456)-w(137|456)+w(137|256)-w(347|256)+w(347|126)-w(457|126)+w(457|123)
    dTemp = w654321 - w731654 + w731652 - w743652 + w743621 - w754621 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //25:w(123|456)-w(237|456)+w(237|145)-w(267|145)+w(267|134)-w(567|134)+w(567|123)
    dTemp = w654321 - w732654 + w732541 - w762541 + w762431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //26:w(123|456)-w(237|456)+w(237|145)-w(267|145)+w(267|135)-w(467|135)+w(467|123)
    dTemp = w654321 - w732654 + w732541 - w762541 + w762531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //27:w(123|456)-w(237|456)+w(237|145)-w(367|145)+w(367|124)-w(567|124)+w(567|123)
    dTemp = w654321 - w732654 + w732541 - w763541 + w763421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //28:w(123|456)-w(237|456)+w(237|145)-w(367|145)+w(367|125)-w(467|125)+w(467|123)
    dTemp = w654321 - w732654 + w732541 - w763541 + w763521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //29:w(123|456)-w(237|456)+w(237|146)-w(257|146)+w(257|134)-w(567|134)+w(567|123)
    dTemp = w654321 - w732654 + w732641 - w752641 + w752431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //30:w(123|456)-w(237|456)+w(237|146)-w(257|146)+w(257|136)-w(457|136)+w(457|123)
    dTemp = w654321 - w732654 + w732641 - w752641 + w752631 - w754631 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //31:w(123|456)-w(237|456)+w(237|146)-w(357|146)+w(357|124)-w(567|124)+w(567|123)
    dTemp = w654321 - w732654 + w732641 - w753641 + w753421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //32:w(123|456)-w(237|456)+w(237|146)-w(357|146)+w(357|126)-w(457|126)+w(457|123)
    dTemp = w654321 - w732654 + w732641 - w753641 + w753621 - w754621 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //33:w(123|456)-w(237|456)+w(237|156)-w(247|156)+w(247|135)-w(467|135)+w(467|123)
    dTemp = w654321 - w732654 + w732651 - w742651 + w742531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //34:w(123|456)-w(237|456)+w(237|156)-w(247|156)+w(247|136)-w(457|136)+w(457|123)
    dTemp = w654321 - w732654 + w732651 - w742651 + w742631 - w754631 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //35:w(123|456)-w(237|456)+w(237|156)-w(347|156)+w(347|125)-w(467|125)+w(467|123)
    dTemp = w654321 - w732654 + w732651 - w743651 + w743521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //36:w(123|456)-w(237|456)+w(237|156)-w(347|156)+w(347|126)-w(457|126)+w(457|123)
    dTemp = w654321 - w732654 + w732651 - w743651 + w743621 - w754621 + w754321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //37:w(123|457)-w(126|457)+w(126|345)-w(167|345)+w(167|234)-w(567|234)+w(567|123)
    dTemp = w754321 - w754621 + w621543 - w761543 + w761432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //38:w(123|457)-w(126|457)+w(126|345)-w(167|345)+w(167|235)-w(467|235)+w(467|123)
    dTemp = w754321 - w754621 + w621543 - w761543 + w761532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //39:w(123|457)-w(126|457)+w(126|345)-w(267|345)+w(267|134)-w(567|134)+w(567|123)
    dTemp = w754321 - w754621 + w621543 - w762543 + w762431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //40:w(123|457)-w(126|457)+w(126|345)-w(267|345)+w(267|135)-w(467|135)+w(467|123)
    dTemp = w754321 - w754621 + w621543 - w762543 + w762531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //41:w(123|457)-w(126|457)+w(126|347)-w(156|347)+w(156|234)-w(567|234)+w(567|123)
    dTemp = w754321 - w754621 + w743621 - w743651 + w651432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //42:w(123|457)-w(126|457)+w(126|347)-w(256|347)+w(256|134)-w(567|134)+w(567|123)
    dTemp = w754321 - w754621 + w743621 - w743652 + w652431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //43:w(123|457)-w(126|457)+w(126|357)-w(146|357)+w(146|235)-w(467|235)+w(467|123)
    dTemp = w754321 - w754621 + w753621 - w753641 + w641532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //44:w(123|457)-w(126|457)+w(126|357)-w(246|357)+w(246|135)-w(467|135)+w(467|123)
    dTemp = w754321 - w754621 + w753621 - w753642 + w642531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //45:w(123|457)-w(136|457)+w(136|245)-w(167|245)+w(167|234)-w(567|234)+w(567|123)
    dTemp = w754321 - w754631 + w631542 - w761542 + w761432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //46:w(123|457)-w(136|457)+w(136|245)-w(167|245)+w(167|235)-w(467|235)+w(467|123)
    dTemp = w754321 - w754631 + w631542 - w761542 + w761532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //47:w(123|457)-w(136|457)+w(136|245)-w(367|245)+w(367|124)-w(567|124)+w(567|123)
    dTemp = w754321 - w754631 + w631542 - w763542 + w763421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //48:w(123|457)-w(136|457)+w(136|245)-w(367|245)+w(367|125)-w(467|125)+w(467|123)
    dTemp = w754321 - w754631 + w631542 - w763542 + w763521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //49:w(123|457)-w(136|457)+w(136|247)-w(156|247)+w(156|234)-w(567|234)+w(567|123)
    dTemp = w754321 - w754631 + w742631 - w742651 + w651432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //50:w(123|457)-w(136|457)+w(136|247)-w(356|247)+w(356|124)-w(567|124)+w(567|123)
    dTemp = w754321 - w754631 + w742631 - w742653 + w653421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //51:w(123|457)-w(136|457)+w(136|257)-w(146|257)+w(146|235)-w(467|235)+w(467|123)
    dTemp = w754321 - w754631 + w752631 - w752641 + w641532 - w764532 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //52:w(123|457)-w(136|457)+w(136|257)-w(346|257)+w(346|125)-w(467|125)+w(467|123)
    dTemp = w754321 - w754631 + w752631 - w752643 + w643521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //53:w(123|457)-w(236|457)+w(236|145)-w(267|145)+w(267|134)-w(567|134)+w(567|123)
    dTemp = w754321 - w754632 + w632541 - w762541 + w762431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //54:w(123|457)-w(236|457)+w(236|145)-w(267|145)+w(267|135)-w(467|135)+w(467|123)
    dTemp = w754321 - w754632 + w632541 - w762541 + w762531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //55:w(123|457)-w(236|457)+w(236|145)-w(367|145)+w(367|124)-w(567|124)+w(567|123)
    dTemp = w754321 - w754632 + w632541 - w763541 + w763421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //56:w(123|457)-w(236|457)+w(236|145)-w(367|145)+w(367|125)-w(467|125)+w(467|123)
    dTemp = w754321 - w754632 + w632541 - w763541 + w763521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //57:w(123|457)-w(236|457)+w(236|147)-w(256|147)+w(256|134)-w(567|134)+w(567|123)
    dTemp = w754321 - w754632 + w741632 - w741652 + w652431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //58:w(123|457)-w(236|457)+w(236|147)-w(356|147)+w(356|124)-w(567|124)+w(567|123)
    dTemp = w754321 - w754632 + w741632 - w741653 + w653421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //59:w(123|457)-w(236|457)+w(236|157)-w(246|157)+w(246|135)-w(467|135)+w(467|123)
    dTemp = w754321 - w754632 + w751632 - w751642 + w642531 - w764531 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //60:w(123|457)-w(236|457)+w(236|157)-w(346|157)+w(346|125)-w(467|125)+w(467|123)
    dTemp = w754321 - w754632 + w751632 - w751643 + w643521 - w764521 + w764321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //61:w(123|467)-w(125|467)+w(125|346)-w(157|346)+w(157|234)-w(567|234)+w(567|123)
    dTemp = w764321 - w764521 + w643521 - w751643 + w751432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //62:w(123|467)-w(125|467)+w(125|346)-w(257|346)+w(257|134)-w(567|134)+w(567|123)
    dTemp = w764321 - w764521 + w643521 - w752643 + w752431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //63:w(123|467)-w(125|467)+w(125|347)-w(156|347)+w(156|234)-w(567|234)+w(567|123)
    dTemp = w764321 - w764521 + w743521 - w743651 + w651432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //64:w(123|467)-w(125|467)+w(125|347)-w(256|347)+w(256|134)-w(567|134)+w(567|123)
    dTemp = w764321 - w764521 + w743521 - w743652 + w652431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //65:w(123|467)-w(135|467)+w(135|246)-w(157|246)+w(157|234)-w(567|234)+w(567|123)
    dTemp = w764321 - w764531 + w642531 - w751642 + w751432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //66:w(123|467)-w(135|467)+w(135|246)-w(357|246)+w(357|124)-w(567|124)+w(567|123)
    dTemp = w764321 - w764531 + w642531 - w753642 + w753421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //67:w(123|467)-w(135|467)+w(135|247)-w(156|247)+w(156|234)-w(567|234)+w(567|123)
    dTemp = w764321 - w764531 + w742531 - w742651 + w651432 - w765432 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //68:w(123|467)-w(135|467)+w(135|247)-w(356|247)+w(356|124)-w(567|124)+w(567|123)
    dTemp = w764321 - w764531 + w742531 - w742653 + w653421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //69:w(123|467)-w(235|467)+w(235|146)-w(257|146)+w(257|134)-w(567|134)+w(567|123)
    dTemp = w764321 - w764532 + w641532 - w752641 + w752431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //70:w(123|467)-w(235|467)+w(235|146)-w(357|146)+w(357|124)-w(567|124)+w(567|123)
    dTemp = w764321 - w764532 + w641532 - w753641 + w753421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //71:w(123|467)-w(235|467)+w(235|147)-w(256|147)+w(256|134)-w(567|134)+w(567|123)
    dTemp = w764321 - w764532 + w741532 - w741652 + w652431 - w765431 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    //72:w(123|467)-w(235|467)+w(235|147)-w(356|147)+w(356|124)-w(567|124)+w(567|123)
    dTemp = w764321 - w764532 + w741532 - w741653 + w653421 - w765421 + w765321;
    if(dTemp <= 0) return 0;
    if(dWei > dTemp) dWei = dTemp;
    
    
    return 0.5*dWei;
}

// EOF
