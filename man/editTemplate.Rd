\name{editTemplate}
\alias{editTemplate}
\title{Edit 3D template}
\usage{
editTemplate(template, fixed, n)
}
\arguments{
  \item{template}{Matrix of template 3D coordinates.}

  \item{fixed}{Number of "fixed" landmark points (non
  surface sliding points)}

  \item{n}{Number of points to be removed}
}
\value{
Function returns a matrix containing the x,y,z coordinates
of the new template landmarks. Function also writes to
working directory 'template.txt' containing the x,y,z
coordinates of the template
}
\description{
An interactive function to remove landmarks from a 3D
template file.
}
\details{
Function edits a 'template.txt' file made by
\code{\link{buildtemplate}}, which must be in current
working directory.  Function overwrites 'template.txt' in
working directory with edited version. Use
read.table("template.txt", header = T).
\subsection{Selection}{ Choosing which landmarks will be
deleted involves landmark selection using a mouse in the
rgl plot window. With a standard 3-button (PC)
buildtemplate uses: \enumerate{ \item the RIGHT mouse
button (primary) to choose points to be deleted (click-drag
a box around landmark), \item the LEFT mouse button
(secondary) is used to rotate mesh, \item the mouse
SCROLLER (third/middle) is used to zoom in and out. } NOTE:
Digitizing functions on MACINTOSH computers using a
standard 3-button mice works as specified. Macs using
platform specific single button mice: \enumerate{ \item
press button to rotate 3D mesh, \item press button while
pressing CONTROL key to choose points to be deleted
(click-drag a box around landmark), \item press button
while pressing OPTION key to zoom in and out.  } Mac mouse
settings need adjustment to acquire zooming functions of
the "middle/third" mouse button/scroller. Some instructions
can be found at \url{http://blog.boastr.net/}. In addition,
freeware providing "Middle Click" functionality is
available at \url{http://magicprefs.com/} for "magicmice"
now standard on many Macintosh machines. }
}
\author{
Erik Otarola-Castillo & Emma Sherratt
}
\keyword{editTemplate}

