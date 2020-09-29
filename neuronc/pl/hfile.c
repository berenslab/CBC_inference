/* Disk file I/O drivers */

/*	Latest Mod	31-Jan-82 		R.G.Smith */

#include <stdio.h>
#include <errno.h>
#include "mdef.h"	/* Global Definitions		  */


#define RDWR 2		/* Read-write for file open */


/************************************/


int prvcreate (name,dnam,ddriv,file,ftyp,mode)
      TEXT   *name;		/* User entered file name */
      TEXT   *dnam;		/* Default file name      */
      TEXT   ddriv;		/* Default disk name  (Ignored in Unix) */
      tfil   *file;		/* File control block     */
      int    ftyp;		/* File type; 0 = process; 1 = text */
      int    mode;		/* Read = 0; Write = 1 */

/* Make a new file on disk, or if one with the
same name already exists, open it instead.

    file->fd 	= file descriptor number;
    file->flen	= length of the file in blocks. 
    file->blk	= current block number in file.
    file->indx	= index into block.
    file->buf	= file buffer for text files.
 
*/

 begin
   int i;
 if (! *name) name = dnam;		/* Default if no name */
 if ((file->fd = open (name, mode)) == -1)   /* Try to find file */
  begin
   if (! mode) file->flen = 0;		/* Don't read non-existent file */
   else if ((file->fd = creat (name,0777)) == -1)  /* Otherwise, try to make */
    begin
     fprintf  (stderr,"Can't make file...");
     file->flen = -1;
    end
    else
     begin
      file->flen = 0;
      close (file->fd);
      file->fd = open (name, RDWR);
     end 
  end
 else
  begin
   file->flen = 1;
   file->blk = 0;
   file->indx = BUFSIZE;
   if (ftyp)
    begin
     if (file->flen == 0)
      for (i=0; i<BUFSIZE; i++)
       file->buf [i] = 0;
   /*  else rbuf (file);		/* Read first block of file */
    end
  end
 return (file->flen);
end


/************************************/

prclose (file)
    pfil *file;

/* Close the file associated with FILE. */


begin
 close (file->fd);
end

/************************************/

BOOL pwrite (file,blk,buf,len)
    pfil *file;
    int  blk,len;
    TEXT *buf;

/* Write out 'len' blocks starting at
address 'buf' onto 'file', at block 'blk'.
Returns T if successful; F otherwise. */

begin
   int i;

  lseek (file->fd, (long) (BUFSIZE * (long)(blk)),0);
  if ((i=write (file->fd,buf, BUFSIZE * len)) == -1)
    fprintf (stderr,"Write error %d at block %d\n",errno,blk);
  return (i ? 0 : 1);
end


/*****************************/

BOOL pread (file,blk,buf,len)
    pfil *file;
    int  blk,len;
    TEXT *buf;

/* Read 'len' blocks starting at
address 'buf' from 'file' at
block 'blk'. */

begin
  lseek (file->fd, (long)(BUFSIZE * (long)(blk)), 0);
  return (read (file->fd, buf, BUFSIZE * len));
end

/********************************/
