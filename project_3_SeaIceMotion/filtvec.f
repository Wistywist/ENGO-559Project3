	program main
	dimension vect(5,179,179),vecta(5,32041)
	character*90 filnm1, filnm2
	character*30 name1,name2
	call getarg(1,filnm1)
	call getarg(2,filnm2)
	open(unit=13,file=filnm1,status='old')
	read(13,999) name1,name2
999	format (2a25)
	read(13,*) ix,jy,iix,jjy,amax
	write(*,*)ix,jy,iix,jjy,amax
	read(*,*)cut
	read(*,*)neigh
	read(*,*)pixels
	do 10 j=1,jy
	do 10 i=1,ix
		read(13,*)(vect(k,i,j),k=1,5)
10	continue
	close (13)
	k=0
	do 60 j=1,jy
	do 50 i=1,ix	
		if(vect(5,i,j).le.cut) go to 50
		iz=0
		ja=j-1
		jb=j+1
		ia=i-1
		ib=i+1
		if(jb.gt.jy)jb=j
		if(ja.lt.1)ja=1
		if(ib.gt.ix)ib=j
		if(ia.lt.1)ia=1
		do 40 jj=ja,jb
		do 30 ii=ia,ib
			if((ii.eq.i).and.(jj.eq.j))go to 30
			if(vect(5,ii,jj).le.cut)go to 30
			xdiff=abs(vect(3,i,j)-vect(3,ii,jj))*amax
			ydiff=abs(vect(4,i,j)-vect(4,ii,jj))*amax
			if((xdiff.le.pixels).and.(ydiff.le.pixels))iz=iz+1
30		continue
40		continue
		if(iz.ge.neigh)then
		k=k+1
		vecta(1,k)=vect(1,i,j)
		vecta(2,k)=vect(2,i,j)
		vecta(3,k)=vect(3,i,j)	
		vecta(4,k)=vect(4,i,j)	
		vecta(5,k)=vect(5,i,j)
		endif
50	continue
60	continue
	write(*,*),k
		
	open(unit=14,file=filnm2,status='unknown')
	write(14,*) name1,name2
	write(14,*)k,1,iix,jjy,amax
	do 80 i=1,k
		write(14,99)(vecta(kk,i),kk=1,5)
80	continue
	close(14)
99	format (9f10.2)
	end

