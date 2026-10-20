from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
"""Combine the monolayer top view (ab plane) and side view (cb plane) into Fig. 2.
b is vertical in both; the side view is scaled so the two unit-cell boxes have the same
b height (596 vs 598.5 px) and shifted so the box edges, hence the atomic rows, align."""
from PIL import Image, ImageDraw, ImageFont
import numpy as np
def load(f):
    im=Image.open(f).convert('RGBA'); bg=Image.new('RGBA',im.size,(255,255,255,255))
    return Image.alpha_composite(bg,im).convert('RGB')
A=load(ROOT/'data/images/TaPdTe5_monolayer_ab.png'); B=load(ROOT/'data/images/TaPdTe5_monolayer_cb.png')
topA,botA=18.0,614.0; topB,botB=77.5,676.0
s=(botA-topA)/(botB-topB)
B2=B.resize((round(B.width*s),round(B.height*s)),Image.LANCZOS)
dy=round(topA-topB*s); gap=50; pad=60
top=min(0,dy); H=max(A.height,B2.height+dy)-top
C=Image.new('RGB',(A.width+gap+B2.width,H+pad),'white')
C.paste(A,(0,pad-top)); C.paste(B2,(A.width+gap,pad+dy-top))
d=ImageDraw.Draw(C)
try: f=ImageFont.truetype('/usr/share/fonts/truetype/crosextra/Carlito-Regular.ttf',40)
except: f=None
ly=pad-top-45; d.text((10,ly),'(a)',fill='black',font=f); d.text((A.width+gap+10,ly),'(b)',fill='black',font=f)
a=np.array(C.convert('L'))<245; ys=np.where(a.any(1))[0]; xs=np.where(a.any(0))[0]
C=C.crop((max(xs[0]-15,0),max(ys[0]-15,0),min(xs[-1]+15,C.width),min(ys[-1]+15,C.height)))
C.save(ROOT/'figures/fig2_monolayer_structure.png',dpi=(300,300))
print('scale',round(s,4),'dy',dy,C.size)
