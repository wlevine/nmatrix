File.open("lapacke.c","w") do |file|
  Dir.glob("lapacke/*/*.c") do |file2|
    file.puts "#include \"#{file2}\""
  end
end
