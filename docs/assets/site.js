// Progressive enhancements: chapters and code remain readable without JavaScript.
const toggle = document.querySelector('.menu-toggle');
const navigation = document.querySelector('#chapter-nav');
toggle?.addEventListener('click', () => {
  const open = navigation.classList.toggle('is-open');
  toggle.setAttribute('aria-expanded', String(open));
});
document.addEventListener('keydown', event => {
  if (event.key === 'Escape' && navigation?.classList.contains('is-open')) {
    navigation.classList.remove('is-open');
    toggle.setAttribute('aria-expanded', 'false');
    toggle.focus();
  }
});
for (const pre of document.querySelectorAll('pre')) {
  const code = pre.querySelector('code');
  if (!code || code.classList.contains('language-text')) continue;
  const button = document.createElement('button');
  button.className = 'copy-button';
  button.textContent = 'Copy';
  button.setAttribute('aria-label', 'Copy code block');
  button.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(code.textContent);
      button.textContent = 'Copied';
    } catch {
      const selection = window.getSelection();
      const range = document.createRange();
      range.selectNodeContents(code);
      selection.removeAllRanges(); selection.addRange(range);
      button.textContent = 'Selected';
    }
    setTimeout(() => { button.textContent = 'Copy'; }, 1800);
  });
  pre.appendChild(button);
}
