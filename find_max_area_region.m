%% 找面积最大的联通区域
function max_region = find_max_area_region(binary)
    % Find the largest connected component in a binary image
    [labeled, num_labels] = bwlabel(binary);
    if num_labels <= 1
        max_region = binary;
        return;
    end
    stats = regionprops(labeled, 'Area');
    areas = [stats.Area];
    [~, max_label] = max(areas);
    max_region = ismember(labeled, max_label) * 255;
end