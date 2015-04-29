subroutine fileman(filename,filelen,filenb,act)

implicit none

integer::filelen
character(len=filelen)::filename
integer::act !2 to re-open, 1 to open, 0 to close, -1 to delete
integer::filenb

integer::io=0

!Opening the file
!==================================================
if (act .eq. 1) then

  open(filenb,file=filename,iostat=io)

  if (io .ne. 0) then

    write (6,'(19a,10a)') "ERROR OPENING FILE ",filename
    write (6,'(9a,i1)') "IOSTAT = ",io
    write (6,'(12a)') "EXITING CODE"
    stop

  end if

else if (act .eq. 2) then

  open(filenb,file=filename,position='append',iostat=io)

  if (io .ne. 0) then

    write (6,'(19a,10a)') "ERROR OPENING ALREADY EXISTING FILE ",filename
    write (6,*) "IOSTAT = ",io
    write (6,'(12a)') "EXITING CODE"
    stop

  end if

!Closing the file
!==================================================
else if (act .eq. 0) then

  close(filenb,iostat=io)

  if (io .ne. 0) then

    write (6,'(19a,10a)') "ERROR CLOSING FILE ",filename
    write (6,'(9a,i1)') "IOSTAT = ",io
    write (6,'(12a)') "EXITING CODE"
    stop

  end if

else if (act .eq. -1) then

  close(filenb, status='delete',iostat=io)

  if (io .ne. 0) then

    write (6,'(20a,10a)') "ERROR DELETING FILE ",filename
    write (6,'(9a,i1)') "IOSTAT = ",io
    write (6,'(12a)') "EXITING CODE"
    stop

  end if

else

  write (6,*) "BAD USE OF fileman, act must be 1, 0 or -1"
  write (6,*) "EXITING CODE"
  stop

endif

end subroutine
