CONFIG = Scheduling_exp_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

tol = 0;

violation_list = false(size(wc_samples, 1), 1);

for exam_id = 1:task_num
    violation_list = violation_list | ((wc_samples(:, exam_id) > st1_deci(exam_id) + tol) & any(wc_samples(:, exam_id + 1:end) ...
        <= st1_deci(exam_id + 1:end)' - tol, 2));
end

fprintf('%d violations found with tolerance %e.\n', sum(violation_list), tol);

tol = 1e-5;

violation_list = false(size(wc_samples, 1), 1);

for exam_id = 1:task_num
    violation_list = violation_list | ((wc_samples(:, exam_id) > st1_deci(exam_id) + tol) & any(wc_samples(:, exam_id + 1:end) ...
        <= st1_deci(exam_id + 1:end)' - tol, 2));
end

fprintf('%d violations found with tolerance %e.\n', sum(violation_list), tol);

tol = 1e-4;

violation_list = false(size(wc_samples, 1), 1);

for exam_id = 1:task_num
    violation_list = violation_list | ((wc_samples(:, exam_id) > st1_deci(exam_id) + tol) & any(wc_samples(:, exam_id + 1:end) ...
        <= st1_deci(exam_id + 1:end)' - tol, 2));
end

fprintf('%d violations found with tolerance %e.\n', sum(violation_list), tol);