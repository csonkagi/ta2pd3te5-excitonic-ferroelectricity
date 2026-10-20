from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
"""Combine the VESTA a- and b-projections of bulk Ta2Pd3Te5 into Fig. 1, panels (a) and (b).
The right panel is scaled so both unit-cell boxes have the same c height (856 vs 777 px)
and shifted so the box edges, hence the atomic rows, align. A red line marks x = a/2,
the mid-plane between the two sheets. Panel labels match Fig. 2."""
from PIL import Image, ImageDraw, ImageFont
import numpy as np
L=Image.open(ROOT/'data/images/TaPdTe5_abulkstructure.png').convert('RGB')
R=Image.open(ROOT/'data/images/TaPdTe5_bbulkstructure.png').convert('RGB')
s=(1336.5-480.5)/(1337.5-560.5)
R2=R.resize((round(R.width*s),round(R.height*s)),Image.LANCZOS)
dy=round(480.5-560.5*s); gap=40
C=Image.new('RGB',(L.width+gap+R2.width,L.height),'white')
C.paste(L,(0,0)); C.paste(R2,(L.width+gap,dy))
a=np.array(C.convert('L'))<245; ys=np.where(a.any(1))[0]; xs=np.where(a.any(0))[0]
x0,y0=max(xs[0]-20,0),max(ys[0]-20,0)
C=C.crop((x0,y0,min(xs[-1]+20,C.width),min(ys[-1]+20,C.height)))
d=ImageDraw.Draw(C); xr=(1037+1680)/2                 # x = a/2 of the right cell box (cropped coords)
d.line([(xr,10),(xr,C.height-10)],fill=(220,0,0),width=5)
# panel labels in a top margin, as in Fig. 2 (font scaled to the image width ratio of the two figures)
top=70; D=Image.new('RGB',(C.width,C.height+top),'white'); D.paste(C,(0,top))
f=ImageFont.truetype('/usr/share/fonts/truetype/crosextra/Carlito-Regular.ttf',round(40*C.width/969))
dd=ImageDraw.Draw(D)
dd.text((10,0),'(a)',fill='black',font=f); dd.text((L.width+gap-x0+10,0),'(b)',fill='black',font=f)
D.save(ROOT/'figures/fig1_bulk_structure.png',dpi=(300,300)); print(D.size)
