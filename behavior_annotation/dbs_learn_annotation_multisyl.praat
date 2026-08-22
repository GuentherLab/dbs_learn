### Score dbs-learn project multisyl accuracy and speech epochs

clearinfo

## GUI (Subject, file type, tiers)
form Select subject, file type, and tiers
    integer starting_file_index 1
    sentence file_name_or_initial_substring trial
    sentence file_extension wav
endform

## Folder Selection
wd$ = chooseDirectory$: "Select the folder containing your sound files"
if wd$ = ""
    exitScript: "No folder selected."
endif
if right$(wd$, 1) <> "/"
    wd$ = wd$ + "/"
endif
outDir$ = wd$

## Normalize inputs
if starting_file_index < 1
    starting_file_index = 1
endif
if file_extension$ = ""
    file_extension$ = "wav"
endif

## Enforced number of tiers
expectedTiers$ = "accuracy transcript disfluency comments unusable_trial difficult_to_score speech_epoch nontarget_sounds_epoch"
logFile$ = wd$ + "accuracy-log.txt"

## Build file list
pattern$ = wd$ + file_name_or_initial_substring$ + "*." + file_extension$
strings = Create Strings as file list: "fileList", pattern$
select Strings fileList
numFiles = Get number of strings
if numFiles = 0
    exitScript: "No files found: " + pattern$
endif
 
## Extract trial numbers and store filenames
for i to numFiles
    filename'i'$ = Get string... i
    ## Extract the number after "trial-"
    name$ = filename'i'$
    hyphenPos = index(name$, "-")
    afterHyphen$ = mid$(name$, hyphenPos + 1, length(name$))
    nextHyphenPos = index(afterHyphen$, "_")
    trialNum'i' = number(left$(afterHyphen$, nextHyphenPos - 1))
endfor
 
## Bubble sort indices by trial number
for i to numFiles
    sortIndex'i' = i
endfor
 
for i to numFiles - 1
    for j to numFiles - i
        a = sortIndex'j'
        jplus = j + 1
        b = sortIndex'jplus'
        if trialNum'a' > trialNum'b'
            sortIndex'j' = b
            sortIndex'jplus' = a
        endif
    endfor
endfor
 
## Build sorted filename array
for i to numFiles
    idx = sortIndex'i'
    sortedFilename'i'$ = filename'idx'$
endfor
 
## Main Loop
ifile = starting_file_index

while ifile >= 1 and ifile <= numFiles
    filename$ = sortedFilename'ifile'$
    Read from file... 'wd$''filename$'
    soundname$ = selected$ ("Sound", 1)

    # Look for TextGrid or enforce number of tiers
    full$ = "'wd$''soundname$'.TextGrid"
    recreate = 0
    if fileReadable (full$)
        Read from file... 'full$'
        select TextGrid 'soundname$'
        nTiers = Get number of tiers
        if nTiers <> 8
            recreate = 1
        else
            # Define expected names
            name1$ = "accuracy"
            name2$ = "transcript"
            name3$ = "disfluency"
            name4$ = "comments"
            name5$ = "unusable_trial"
            name6$ = "difficult_to_score"
            name7$ = "speech_epoch"
            name8$ = "nontarget_sounds_epoch"

            # Query current tier names as variables
            tier1$ = Get tier name... 1
            tier2$ = Get tier name... 2
            tier3$ = Get tier name... 3
            tier4$ = Get tier name... 4
            tier5$ = Get tier name... 5
            tier6$ = Get tier name... 6
            tier7$ = Get tier name... 7
            tier8$ = Get tier name... 8

            # Now, do the comparison:
            if tier1$ <> name1$ or tier2$ <> name2$ or tier3$ <> name3$ or tier4$ <> name4$ or tier5$ <> name5$ or tier6$ <> name6$ or tier7$ <> name7$ or tier8$ <> name8$
                recreate = 1
            endif
        endif
        if recreate = 1
            Remove
            select Sound 'soundname$'
            To TextGrid... "'expectedTiers$'" ""
        endif
    else
        select Sound 'soundname$'
        To TextGrid... "'expectedTiers$'" ""
    endif

    # Open editor
    select Sound 'soundname$'
    plus TextGrid 'soundname$'
    View & Edit

    # Pause for scoring, with navigation
    beginPause: "Scoring file 'ifile' of 'numFiles': 'filename$'"
        comment: "Current file: 'filename$'"
        optionMenu: "Navigation", 1
            option: "Continue to next file"
            option: "Jump to specific file"
            option: "Skip ahead"
        optionMenu: "Select file", 1
            for i to numFiles
                option: sortedFilename'i'$
            endfor
        optionMenu: "Skip amount", 1
            option: "1"
            option: "10"
            option: "25"
    clicked = endPause: "Continue", "Quit", 1

    # Save TextGrid
    select TextGrid 'soundname$'
    Save as text file: wd$ + soundname$ + ".TextGrid"

    # Extract labels
    accuracy$               = Get label of interval... 1 1
    transcript$             = Get label of interval... 2 1
    disfluency$             = Get label of interval... 3 1
    comments$               = Get label of interval... 4 1
    unusable_trial$         = Get label of interval... 5 1
    difficult_to_score$     = Get label of interval... 6 1
    speech_epoch$           = Get label of interval... 7 1
    nontarget_sounds_epoch$ = Get label of interval... 8 1

    # Append log line
    fileappend 'logFile$' 'filename$' \t 'accuracy$' \t 'transcript$' \t 'disfluency$' \t 'comments$' \t 'unusable_trial$' \t 'difficult_to_score$' \t 'speech_epoch$' \t 'nontarget_sounds_epoch$' \n

    # Cleanup objects
    select TextGrid 'soundname$'
    plus Sound 'soundname$'
    Remove
    clearinfo
    select Strings fileList

    # Quit
    if clicked = 2
        select Strings fileList
        Remove
        exitScript: "Stopped by user."
    endif

    # Navigation logic
    if navigation = 1
        ifile = ifile + 1
    elsif navigation = 2
        ifile = select_file
    elsif navigation = 3
        if skip_amount = 1
            skip = 1
        elsif skip_amount = 2
            skip = 10
        else
            skip = 25
        endif
        ifile = ifile + skip
    endif
endwhile

## Cleanup
select Strings fileList
Remove
printline All done! Scored 'numFiles' files. Log saved to 'logFile$'.