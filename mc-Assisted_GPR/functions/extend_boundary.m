function farVal = extend_boundary(maxOrMinVal, multiplier)
    farVal = sign(maxOrMinVal) * abs(maxOrMinVal) * multiplier;
end

