%==============================================================================
% THESIS CLASS OPTIONS
%==============================================================================
% Define as class file
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{mathnotes}[1995/12/01 v1.0 Math Notes class]

% Provide warning if user specifies an option since none have been defined
\DeclareOption*{\typeout{Unknown option \CurrentOption.}}

% Pass options to the report class before loading
\PassOptionsToClass{final}{report}
\PassOptionsToClass{letterpaper}{report}
\PassOptionsToClass{12pt}{report}

% Process class options
\ProcessOptions

% Load 'report' class to act as starting point
\LoadClass{article}

%==============================================================================
% MARGINS
%==============================================================================
\setlength{\oddsidemargin}{0in}
\setlength{\evensidemargin}{0in}
\setlength{\topmargin}{0in}
\setlength{\headsep}{0in}
\setlength{\headheight}{0in}
\setlength{\textheight}{9in}
\setlength{\textwidth}{6.5in}

%==============================================================================
% THEOREM ENVIRONMENTS
%==============================================================================
\RequirePackage[most]{tcolorbox}
%\RequirePackage{chngcntr}
%\RequirePackage{amsthm}

% counters
%\newcounter{theorem}
%\newcounter{lemma}
%\newcounter{corollary}
%\counterwithin{theorem}{section}
%\counterwithin{lemma}{section}
%\counterwithin{corollary}{section}

\makeatletter

% environment for theorems
\newtcolorbox[auto counter,number within=section]{theorem}[2][]{
breakable,
enhanced,
colback=blue!05,
colframe=blue,
top=\baselineskip,
enlarge top by=\topsep,
overlay unbroken and first={
  \node[xshift=10pt,thick,draw=blue,fill=blue!15,rounded corners,anchor=west] at (frame.north west) %
  {\strut{\bfseries Theorem~\thetcbcounter}\if#2\@empty\relax\relax\else~(#2)\fi};
  },
#1
}

% environment for lemmas
\newtcolorbox[auto counter,number within=section]{lemma}[2][]{
breakable,
enhanced,
colback=green!05,
colframe=green!90!black,
top=\baselineskip,
enlarge top by=\topsep,
overlay unbroken and first={
  \node[xshift=10pt,thick,draw=green!90!black,fill=green!20,rounded corners,anchor=west] at (frame.north west) %
  {\strut{\bfseries Lemma~\thetcbcounter}\if#2\@empty\relax\relax\else~(#2)\fi};
  },
#1
}

% environment for Corollaries
\newtcolorbox[auto counter,number within=section]{corollary}[2][]{
breakable,
enhanced,
colback=red!05,
colframe=red!90!black,
top=\baselineskip,
enlarge top by=\topsep,
overlay unbroken and first={
  \node[xshift=10pt,thick,draw=red!90!black,fill=red!20,rounded corners,anchor=west] at (frame.north west) %
  {\strut{\bfseries Corollary~\thetcbcounter}\if#2\@empty\relax\relax\else~(#2)\fi};
  },
#1
}

% environment for Definitions
\newtcolorbox[auto counter,number within=section]{definition}[2][]{
breakable,
enhanced,
colback=orange!05,
colframe=orange!90!black,
top=\baselineskip,
enlarge top by=\topsep,
overlay unbroken and first={
  \node[xshift=10pt,thick,draw=orange!90!black,fill=orange!20,rounded corners,anchor=west] at (frame.north west) %
  {\strut{\bfseries Definition~\thetcbcounter}\if#2\@empty\relax\relax\else~(#2)\fi};
  },
#1
}

% environment for Propositions
\newtcolorbox[auto counter,number within=section]{proposition}[2][]{
breakable,
enhanced,
colback=black!05,
colframe=black!90,
top=\baselineskip,
enlarge top by=\topsep,
overlay unbroken and first={
  \node[xshift=10pt,thick,draw=black!90,fill=black!20,rounded corners,anchor=west] at (frame.north west) %
  {\strut{\bfseries Proposition~\thetcbcounter}\if#2\@empty\relax\relax\else~(#2)\fi};
  },
#1
}

\makeatother

% other math environments
\newenvironment{proof}[1][Proof]{\begin{trivlist}
\item[\hskip \labelsep {\bfseries #1}]}{\end{trivlist}}
\newenvironment{example}[1][Example]{\begin{trivlist}
\item[\hskip \labelsep {\bfseries #1}]}{\end{trivlist}}
\newenvironment{remark}[1][Remark]{\begin{trivlist}
\item[\hskip \labelsep {\bfseries #1}]}{\end{trivlist}}

% QED symbol command
\newcommand{\qed}{\nobreak \ifvmode \relax \else
      \ifdim\lastskip<1.5em \hskip-\lastskip
      \hskip1.5em plus0em minus0.5em \fi \nobreak
      \vrule height0.75em width0.5em depth0.25em\fi}