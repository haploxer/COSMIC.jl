macro test(file)
    if test
        return esc(symbol("test_", file))
    else
        return esc(file)
    end
end
