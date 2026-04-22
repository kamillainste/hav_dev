# HAV Development Project

A collaborative project for viral genome analysis using NextClade.

## 🚀 Getting Started with GitHub Codespaces

This project is configured to work seamlessly in GitHub Codespaces, with NextClade pre-installed.

### For Project Members

1. **Open in Codespaces:**
   - Go to https://github.com/kamillainste/hav_dev
   - Click the green `Code` button
   - Select the `Codespaces` tab
   - Click `Create codespace on main` (or select existing codespace)

2. **Wait for Setup:**
   - The environment will automatically install NextClade and dependencies
   - This takes 2-3 minutes on first launch

3. **Start Working:**
   - NextClade is ready to use: `nextclade --help`
   - All changes are automatically synced via Git

### Using NextClade

```bash
# Check version
nextclade --version

# Run analysis (example)
nextclade run --input-dataset <dataset> <input.fasta>

# Get help
nextclade --help
```

## 📁 Project Structure

- `script.R` - R analysis scripts
- `.devcontainer/` - Codespaces configuration (auto-setup)

## 🤝 Collaboration

- All team members can create their own Codespace from this repository
- Changes are synced via Git commits and pushes
- No local installation required - everything runs in the cloud

## 📝 Notes

This project uses GitHub Codespaces to bypass local installation restrictions and enable seamless collaboration.
