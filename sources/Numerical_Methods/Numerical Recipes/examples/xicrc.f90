	PROGRAM xicrc
	USE nrtype
	USE nr
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(256) :: lin
	INTEGER(I2B) :: i1,i2,j
	INTEGER(I4B) :: n
	do
		write(*,*) 'ENTER LENGTH,STRING: '
		read(*,'(i3,80a1)',END=99) n,(lin(j),j=1,n)
		write(*,*) (lin(j),j=1,n)
		i1=icrc(0_I2B,lin(1:n),0_I2B,1)
		lin(n+1)=char(iand(255_I2B,ishft(i1,-8)))
		lin(n+2)=char(iand(i1,255_I2B))
		i2=icrc(i1,lin(1:n+2),0_I2B,1)
		write(*,'(''    XMODEM: String CRC, Packet CRC= '',z4,2x,z4)')&
			i1,i2
		i1=icrc(i2,lin(1:n),255_I2B,-1)
		lin(n+1)=char(iand(255_I2B,not(iand(i1,255_I2B))))
		lin(n+2)=char(iand(255_I2B,not(ishft(i1,-8))))
		i2=icrc(i1,lin(1:n+2),255_I2B,-1)
		write(*,'(''      X.25: String CRC, Packet CRC= '',z4,2x,z4)')&
			i1,i2
		i1=icrc(0_I2B,lin(1:n),0_I2B,-1)
		lin(n+1)=char(iand(i1,255_I2B))
		lin(n+2)=char(iand(255_I2B,ishft(i1,-8)))
		i2=icrc(i1,lin(1:n+2),0_I2B,-1)
		write(*,'('' CRC-CCITT: String CRC, Packet CRC= '',z4,2x,z4)')&
			i1,i2
	end do
99	write (*,*) 'Normal completion'
	END PROGRAM xicrc
