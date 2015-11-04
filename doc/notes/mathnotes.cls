%==============================================================================
% THESIS CLASS OPTIONS
%==============================================================================
% Define as class file
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{mathnotes}[1995/12/01 v1.0 Math Notes class]

% Provide warning if user specifies an option since none have been defined
\DeclareOption*{\typeout{Unknown option \CurrentOption.}}

% Pass options to the report class before loading
\PassOptionsToClass{final}{article}
\PassOptionsToClass{letterpaper}{article}
\PassOptionsToClass{10pt}{article}

% Process class options
\ProcessOptions

% Load class to act as starting point
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

