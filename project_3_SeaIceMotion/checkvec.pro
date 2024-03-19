PRO checkvec,vecfile
;device,decompose=0
!p.background=255b
!p.color=0
!order=1

xsz=960
ysz=720

x11=0
x22=xsz-1
y11=0
y22=ysz-1
ddx=x22-x11+1
ddy=y22-y11+1
;window,0,xsize=xsz,ysize=ysz
set_plot,'ps'


mapfile='barrow.map'
map=bytarr(xsz,ysz)
openr,1,mapfile
readu,1,map
close,1
map=255b-map

print,vecfile


tv,map
plot,[x11,x22],[y11,y22],/nodata,/xst,/yst, $
    xmargin=[0,0],ymargin=[0,0],/noerase, $
    xticks=1,yticks=1,xrange=[x11,x22], $
    yrange=[ysz-y22,ysz-y11],ticklen=0.,charsize=.01

dum=''
inp=fltarr(5)
openr,1,vecfile
readf,1,dum
readf,1,inp
vect=fltarr(5,inp(0),inp(1))
readf,1,vect
close,1
vectu=reform(vect(2,*,*))
vectv=reform(vect(3,*,*))
vectx=reform(vect(0,*,*))
vecty=reform(vect(1,*,*))
vectc=reform(vect(4,*,*))
 
cut=.0001
maxvec=10
r = .3				;len of arrow head
length=1.0
angle = 22.5 * !dtor		;Angle of arrowhead
st = r * sin(angle)		;sin 22.5 degs * length of head
ct = r * cos(angle)
for j=0,inp(1)-1  do for i=0,inp(0)-1 do if vectc(i,j) gt cut then begin
	x0 = vectx(i,j)
	dx = vectu(i,j)*length
	x1 = (x0 + dx)
	y0 = (ysz-vecty(i,j))
	dy = vectv(i,j)*length
	y1 = (y0 + dy)
	plots,[x0,x1,x1-(ct*dx-st*dy),x1,x1-(ct*dx+st*dy)], $
	      [y0,y1,y1-(ct*dy+st*dx),y1,y1-(ct*dy-st*dx)],color=0,noclip=0
endif

x0 = 22
dx = .5*maxvec*length
x1 = (x0 + dx)
y0 = (25)
dy = 0.0*length
y1 = (y0 + dy)
plots,[x0,x1,x1-(ct*dx-st*dy),x1,x1-(ct*dx+st*dy)],color=0, $
      [y0,y1,y1-(ct*dy+st*dx),y1,y1-(ct*dy-st*dx)]
xyouts,x0,y0-15,strcompress(maxvec/2,/remove_all)+' cm/sec',/data,color=0,charsize=1.2
x0 = 85
dx = maxvec*length
x1 = (x0 + dx)
y0 = (25)
dy = 0.0*length
y1 = (y0 + dy)
plots,[x0,x1,x1-(ct*dx-st*dy),x1,x1-(ct*dx+st*dy)],color=0, $
      [y0,y1,y1-(ct*dy+st*dx),y1,y1-(ct*dy-st*dx)]
xyouts,x0,y0-15,strcompress(maxvec,/remove_all)+' cm/sec',/data,color=0,charsize=1.2
outf=strmid(vecfile,34,10)
yr=strmid(outf,0,4)
mon=strmid(outf,4,3)
xyouts,.01,.91,mon+'  '+yr,/normal,charsize=1.5

end
