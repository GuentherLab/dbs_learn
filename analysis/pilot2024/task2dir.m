%%%% translate task name used in matlab variables (no dashes) to task name used in early versions of filenames (some include dashes)

function dir_string = task2dir(task_string)

switch task_string
    case  'assess',
        dir_string = 'assess';
    case 'pretest',
        dir_string = 'pretest';
    case 'trainA',
        dir_string = 'train-a';
    case'trainB',
        dir_string = 'train-b';
    case'test'
        dir_string = 'test';
    otherwise
        error(['task name ' dir_string ' not recognized'])
end
